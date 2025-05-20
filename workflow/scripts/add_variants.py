import os
import sys
import datetime
import pandas as pd
import numpy as np
from itertools import chain, combinations
import concurrent.futures
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import IntervalTree  # pour l'arbre d'intervalles

# ------------------------------------------------------------
# 1) Configuration du log : redirige stdout & stderr vers un fichier
# ------------------------------------------------------------
workdir = os.getcwd()
log_dir = os.path.join(workdir, "logs")
os.makedirs(log_dir, exist_ok=True)
log_path = os.path.join(log_dir, "apply_variants.log")

# Ouvre le fichier de log en mode line‑buffered
log_fh = open(log_path, "w", buffering=1)

def log(msg: str):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_fh.write(f"{ts} - {msg}\n")
    log_fh.flush()

# Redirige stdout & stderr vers le log
sys.stdout = log_fh
sys.stderr = log_fh

# ------------------------------------------------------------
# 2) Fonctions du pipeline
# ------------------------------------------------------------

def read_vcf(vcf_file: str) -> pd.DataFrame:
    """
    Lit un VCF simple (pas bgzippé) et retourne un DataFrame avec
    colonnes ['chromosome','position','ref_nucleotide','alt_nucleotide']
    """
    mutations = []
    with open(vcf_file, 'r') as v:
        for line in v:
            if line.startswith('#'):
                continue
            chrom, pos, _, ref, alt, *rest = line.strip().split('\t')
            mutations.append((chrom, int(pos), ref, alt))
    df = pd.DataFrame(mutations,
                      columns=['chromosome','position','ref_nucleotide','alt_nucleotide'])
    df['position'] = df['position'].astype(np.uint32)
    log(f"1) VCF chargé avec succès → {len(df)} mutations")
    return df

# Génère toutes les combinaisons non vides
def all_subsets(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

# Variables globales partagées par les workers
global_genome = None
global_gtf_exons = None
global_trees = None

def init_worker(fasta_path: str, exon_parquet_path: str):
    """Initialise le worker : charge le génome, le DataFrame d'exons et construit les IntervalTrees"""
    global global_genome, global_gtf_exons, global_trees
    # 1) Genome en mémoire
    global_genome = { rec.id.split('.')[0]: rec
                      for rec in SeqIO.parse(fasta_path, "fasta") }
    # 2) DataFrame d'exons
    global_gtf_exons = pd.read_parquet(exon_parquet_path)
    # 3) Construire un IntervalTree par chromosome
    global_trees = {}
    for exon in global_gtf_exons.itertuples():
        chrom = exon.chromosome
        start = exon.start
        end = exon.end
        if start > end:
            start, end = end, start  # inverse si start > end
        end += 1  # intervaltree half-open
        tree = global_trees.setdefault(chrom, IntervalTree())
        tree.addi(start, end, exon.transcript_id)
    log(f"Clés global_trees : {list(global_trees.keys())}")
    chrom = "1"
    df1 = global_gtf_exons[global_gtf_exons["chromosome"] == chrom]
    log(f"Exons chr{chrom}: start min={df1['start'].min()}, start max={df1['start'].max()}, end min={df1['end'].min()}, end max={df1['end'].max()}")

# Récupère et assemble les exons d'un transcript
def _print_transcript_exons(transcript_id):
    global global_genome, global_gtf_exons
    exons = global_gtf_exons[global_gtf_exons['transcript_id'] == transcript_id]
    if exons.empty or transcript_id not in global_genome:
        return None, None
    strand = exons.iloc[0]['strand']
    sorted_exons = exons.sort_values('start', ascending=(strand == '+'))
    seq = str(global_genome[transcript_id].seq)
    mapping, recon = {}, []
    offset = 0
    for _, exon in sorted_exons.iterrows():
        length = exon['end'] - exon['start'] + 1
        coords = (range(exon['start'], exon['end']+1)
                  if strand=='+' else
                  range(exon['end'], exon['start']-1, -1))
        for i, gpos in enumerate(coords):
            mapping[gpos] = offset + i
        recon.append(seq[offset:offset+length])
        offset += length
    return "".join(recon), mapping

# Applique les mutations sur la séquence du transcript
def _apply_mutations_to_transcript(transcript_id, seq, mapping, muts):
    muts_in = [m for m in muts if m['position'] in mapping]
    applied_mutations = []
    for m in muts_in:
        m['transcript_pos'] = mapping[m['position']]
    muts_in.sort(key=lambda x: x['transcript_pos'])
    new_seq, shift = seq, 0
    for m in muts_in:
        p   = m['transcript_pos'] + shift
        ref = m['ref_nucleotide']
        alt = m['alt_nucleotide']
        if new_seq[p:p+len(ref)] != ref:
            log(f"[WARN] Mutation ignorée pour {transcript_id} à pos transcript {p+1}: attendu {ref}, trouvé {new_seq[p:p+len(ref)]}")
            continue
        new_seq = new_seq[:p] + alt + new_seq[p+len(ref):]
        shift += len(alt) - len(ref)
        applied_mutations.append(m)
    return new_seq, applied_mutations


# Traite un transcript dans un worker
def _process_transcript_worker(transcript_id, mutations):
    global global_genome, global_gtf_exons
    log(f"Début du traitement pour {transcript_id} ({len(mutations)} mutations)")
    exons = global_gtf_exons[global_gtf_exons['transcript_id'] == transcript_id]
    # skip lncRNA
#    if exons.empty or exons.iloc[0]['gene_biotype'] == 'lncRNA':
#        return []
    seq, mapping = _print_transcript_exons(transcript_id)
    if seq is None:
        return []
    out = []
    for combo in all_subsets(mutations):
        combo = list(combo)
        if any(m['position'] not in mapping for m in combo):
            continue
        ms, applied = _apply_mutations_to_transcript(transcript_id, seq, mapping, combo)

        if not applied:
            continue  

        signature = "_".join(f"{m['ref_nucleotide']}({m['alt_nucleotide']})" for m in applied)
        positions = ";".join(str(mapping[m['position']] + 1) for m in applied)


        new_id = f"{transcript_id}_mut_{signature}_pos={positions}"
        new_description = global_genome[transcript_id].description
        log(f"[DEBUG] type(ms)={type(ms)} value(ms)={ms}")
        rec = SeqRecord(ms if isinstance(ms, Seq) else Seq(ms), id=new_id, description=new_description)
        out.append(rec)
    return out 


# Génère les transcrits mutés en parallèle en utilisant l'IntervalTree
def _generate_mutated_transcripts_parallel(
    fasta, vcf_df, exon_parquet, max_workers=12
):
    log("2) Chargement des exons pour le main…")
    # Reconstruit les arbres pour le thread principal
    init_worker(fasta, exon_parquet)
    for chrom, tree in global_trees.items():
        intervals = list(tree)
        if intervals:
            min_start = min(iv.begin for iv in intervals)
            max_end = max(iv.end for iv in intervals)
            log(f"IntervalTree {chrom}: {len(intervals)} intervalles, start min={min_start}, end max={max_end}")
        else:
            log(f"IntervalTree {chrom}: AUCUN intervalle.")

    log(f"   → {len(global_gtf_exons)} exons chargés.")
    # 3) Mapping transcript → mutations via IntervalTree
    tx_muts = {}
    for mut in vcf_df.to_dict('records'):
        chrom = mut['chromosome']
        pos   = mut['position']
        if chrom not in global_trees:
            continue
        for iv in global_trees[chrom][pos]:
            tx = iv.data
            tx_muts.setdefault(tx, []).append(mut)
    total = len(tx_muts)
    log(f"3) Mapping construit → {total} transcripts à traiter.")
    # 4) Exécution parallèle
    log("4) Lancement du ProcessPoolExecutor…")
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=init_worker,
        initargs=(fasta, exon_parquet)
    ) as exe:
        futures = {
            exe.submit(_process_transcript_worker, tx, muts): tx
            for tx, muts in tx_muts.items()
        }
        for i, fut in enumerate(concurrent.futures.as_completed(futures), start=1):
            log(f"[{i}/{total}] Transcript {futures[fut]} traité.")
            for rec in fut.result():
                yield rec
                

# Ecrit les transcrits mutés par batch pour économiser la mémoire
def save_mutated_transcripts_in_batches_parallel(
    fasta, vcf_df, exon_parquet, out_fasta,
    batch_size=10, max_workers=12
):
    batch = []
    with open(out_fasta, "w") as out:
        for rec in _generate_mutated_transcripts_parallel(
            fasta, vcf_df, exon_parquet, max_workers
        ):
            batch.append(rec)
            if len(batch) >= batch_size:
                SeqIO.write(batch, out, "fasta")
                out.flush(); os.fsync(out.fileno())
                batch.clear()
        if batch:
            SeqIO.write(batch, out, "fasta")
            out.flush(); os.fsync(out.fileno())

# ------------------------------------------------------------
# 3) Fonction main()
# ------------------------------------------------------------
def main():
    log("=== DÉMARRAGE DU PIPELINE ===")
    vcf_file          = snakemake.input.vcf
    transcripts_fasta = snakemake.input.transcripts_fasta
    mutated_output_fa = snakemake.output.mutated_output_fasta
    exon_parquet      = snakemake.input.exon_parquet
    
#    vcf_file = "20QC_variant_tronque.vcf"
#    transcripts_fasta = "breast_cancer/pickle/gencode.v47.transcripts.fa"
#    mutated_output_fa  = "mutated_transcripts.fa"
#    exon_parquet = "exon_data.parquet"  
#    print("bob")
    try:
        df = read_vcf(vcf_file)
        save_mutated_transcripts_in_batches_parallel(
            transcripts_fasta, df, exon_parquet,
            mutated_output_fa,
            batch_size=10,
            max_workers=12
        )
        log("=== PIPELINE TERMINÉ AVEC SUCCÈS ===")
    except Exception as e:
        log(f"!!! ERREUR CRITIQUE : {e}")
    finally:
        log_fh.close()

if __name__ == "__main__":
    main()

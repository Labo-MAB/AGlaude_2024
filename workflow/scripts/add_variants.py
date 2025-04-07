import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import logging
from itertools import chain, combinations
import concurrent.futures
import cProfile, pstats


# Configuration du logger
logging.basicConfig(
    filename="logs/apply_variants.log",
    level=logging.DEBUG, # a voir si on remplace DEBUG -> WARNING
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def read_vcf(vcf_file):
    """Charge le fichier VCF en DataFrame en extrayant uniquement les positions et les nucléotides."""
    mutations = []
    try:
        with open(vcf_file, 'r') as vcf:
            for line in vcf:
                if line.startswith('##') or line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                chromosome = fields[0]
                position = int(fields[1])
                ref_nucleotide = fields[3]
                alt_nucleotide = fields[4]
                mutations.append((chromosome, position, ref_nucleotide, alt_nucleotide))
        
        mutations_df = pd.DataFrame(mutations, columns=['chromosome', 'position', 'ref_nucleotide', 'alt_nucleotide'])
        mutations_df['position'] = mutations_df['position'].astype(np.uint32)
        logging.info("VCF chargé avec succès.")
        print("Mutations DataFrame:\n", mutations_df)
        return mutations_df
    except Exception as e:
        logging.error(f"Erreur lors de la lecture du VCF: {e}")
        raise

def load_exons(exon_parquet):
    """Charge les exons depuis un fichier Parquet pour éviter un recalcul."""
    try:
        GTF_exons = pd.read_parquet(exon_parquet)
        GTF_exons['exon_interval'] = pd.IntervalIndex.from_arrays(GTF_exons['start'], GTF_exons['end'], closed="both")
        logging.info("Exons chargés avec succès.")
        return GTF_exons
    except Exception as e:
        logging.error(f"Erreur lors du chargement des exons: {e}")
        print(f"Erreur lors du chargement des exons: {e}", flush=True)
        raise

def _print_transcript_exons(transcript_id, GTF_exons, genome):
    """
    Reconstruit la séquence transcriptomique et établit un mapping entre positions génomiques et transcriptomiques.
    
    Retourne :
      - reconstructed_seq : séquence complète du transcript
      - genomic_to_transcript : dictionnaire {position génomique -> position transcriptomique}
    """
    exons_info = GTF_exons[GTF_exons['transcript_id'] == transcript_id]
    if exons_info.empty:
        msg = f"Aucun exon trouvé pour le transcript {transcript_id}."
        logging.warning(msg)
        print(msg)
        return None, None

    strand = exons_info.iloc[0]['strand']
    sorted_exons = exons_info.sort_values(by='start', ascending=(strand == '+'))
    
    try:
        transcript_seq = str(genome[transcript_id].seq)
    except KeyError as e:
        logging.error(f"Transcript {transcript_id} non trouvé dans le FASTA: {e}")
        return None, None

    genomic_to_transcript = {}
    reconstructed = []
    cumulative = 0  # Position dans le transcript reconstruit

    print(f"\n### DEBUG : Reconstruction du transcript {transcript_id}")
    for i, (_, exon) in enumerate(sorted_exons.iterrows(), start=1):
        exon_len = exon['end'] - exon['start'] + 1
        if strand == '+':
            for j, genomic_pos in enumerate(range(exon['start'], exon['end'] + 1)):
                genomic_to_transcript[genomic_pos] = cumulative + j
        else:
            for j, genomic_pos in enumerate(range(exon['end'], exon['start'] - 1, -1)):
                genomic_to_transcript[genomic_pos] = cumulative + j
        exon_seq = transcript_seq[cumulative: cumulative + exon_len]
        reconstructed.append(exon_seq)
        logging.debug(f"Exon {i} - Genomic: {exon['start']}-{exon['end']}, Transcript: {cumulative}-{cumulative+exon_len-1}")
        print(f"Exon {i} - Coordonnées génomiques: {exon['start']} - {exon['end']}, positions transcriptomiques: {cumulative} - {cumulative+exon_len-1}")
        print(f"Séquence exon: {exon_seq}\n")
        cumulative += exon_len

    reconstructed_seq = "".join(reconstructed)
    logging.info(f"Transcript {transcript_id} reconstruit, longueur {len(reconstructed_seq)}.")
    print("Transcript complet reconstruit :", reconstructed_seq)
    print("Longueur =", len(reconstructed_seq))
    return reconstructed_seq, genomic_to_transcript

def _apply_mutations_to_transcript(transcript_seq, genomic_to_transcript, mutations):
    """
    Applique une liste de mutations à une séquence transcriptomique.
    """
    mutation_list = []
    for m in mutations:
        genomic_pos = m['position']
        if genomic_pos in genomic_to_transcript:
            m["transcript_pos"] = genomic_to_transcript[genomic_pos]
            mutation_list.append(m)
        else:
            msg = f"Position {genomic_pos} non trouvée dans le mapping pour la mutation {m}."
            logging.warning(msg)
            print(msg)
    
    mutation_list.sort(key=lambda x: x["transcript_pos"])
    offset = 0
    new_seq = transcript_seq
    for m in mutation_list:
        pos = m["transcript_pos"] + offset
        ref = m["ref_nucleotide"]
        alt = m["alt_nucleotide"]
        if new_seq[pos: pos + len(ref)] != ref:
            msg = f"**Mismatch** à la position {pos} dans le transcript: attendu {ref}, trouvé {new_seq[pos: pos + len(ref)]}"
            logging.error(msg)
            print(msg)
        new_seq = new_seq[:pos] + alt + new_seq[pos+len(ref):]
        offset += len(alt) - len(ref)
    return new_seq

def all_subsets(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

############################
# Traitement d'un transcript
############################

def _process_transcript_worker(transcript_id, mutations, GTF_exons, genome):
    """
    Traite un transcript en :
      - reconstruisant sa séquence
      - calculant la position transcriptomique pour chaque mutation
      - générant toutes les combinaisons de mutations applicables
      - retournant une liste de SeqRecord mutés
    """
    exons_info = GTF_exons[GTF_exons['transcript_id'] == transcript_id]
    if exons_info.empty:
        logging.warning(f"Aucun exon trouvé pour le transcript {transcript_id}.")
        return []
    if exons_info.iloc[0]['gene_biotype'] == 'lncRNA': # filtred out les lncRNA et transcrit mis fichier a part
        logging.info(f"Transcript {transcript_id} est un lncRNA. Traitement des mutations ignoré.")
        return []
    
    transcript_seq, genomic_to_transcript = _print_transcript_exons(transcript_id, GTF_exons, genome)
    if transcript_seq is None or genomic_to_transcript is None:
        logging.error(f"Impossible de reconstruire le transcript {transcript_id}")
        return []
    
    # Affecter la position dans le transcript à chaque mutation
    for m in mutations:
        genomic_pos = m['position']
        if genomic_pos in genomic_to_transcript:
            m["transcript_pos"] = genomic_to_transcript[genomic_pos]
        else:
            logging.warning(f"Position {genomic_pos} non trouvée pour le transcript {transcript_id}")
            m["transcript_pos"] = None
    
    seq_records = []
    for comb in all_subsets(mutations):
        comb_list = list(comb)
        if any(m["transcript_pos"] is None for m in comb_list):
            continue
        comb_list.sort(key=lambda x: x["transcript_pos"])
        new_seq = _apply_mutations_to_transcript(transcript_seq, genomic_to_transcript, comb_list)
        signature = "_".join([f"{m['ref_nucleotide']}({m['alt_nucleotide']})" for m in comb_list])
        # Création du SeqRecord muté
        if transcript_id in genome:
            original_record = genome[transcript_id]
            new_id = f"{transcript_id}_mut_({signature})"
            new_description = original_record.description
        else:
            new_id = f"{transcript_id}_mut_({signature})"
            new_description = f"Mutated transcript {transcript_id}"
        seq_record = SeqRecord(Seq(new_seq), id=new_id, description=new_description)
        seq_records.append(seq_record)
        logging.info(f"Transcript {transcript_id} - combinaison {signature} appliquée.")
        print(f"Transcript {transcript_id} - combinaison {signature} appliquée.")
    return seq_records

##########################################
# Générateur parallèle de transcrits mutés
##########################################

def _generate_mutated_transcripts_parallel(transcripts_fasta, mutations_df, GTF_exons, max_workers=32):
    genome = {
        record.id.split('.')[0]: record
        for record in SeqIO.parse(transcripts_fasta, "fasta")
    }
    # Association des mutations aux transcrits via les exons
    transcript_mutations = {}
    for idx, mutation in mutations_df.iterrows():
        chromosome = mutation['chromosome']
        position = mutation['position']
        matching_exons = GTF_exons[
            (GTF_exons['chromosome'] == chromosome) &
            (GTF_exons['start'] <= position) &
            (GTF_exons['end'] >= position)
        ]
        if matching_exons.empty:
            logging.warning(f"Aucun exon ne correspond à la mutation à la position {position}.")
            continue
        for _, exon in matching_exons.iterrows():
            transcript_id = exon['transcript_id']
            if transcript_id not in genome:
                logging.error(f"Transcript {transcript_id} non trouvé dans le FASTA.")
                continue
            transcript_mutations.setdefault(transcript_id, []).append(mutation.to_dict())
    
    # Traiter chaque transcript en parallèle
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_process_transcript_worker, tid, muts, GTF_exons, genome): tid
            for tid, muts in transcript_mutations.items()
        }
        for future in concurrent.futures.as_completed(futures):
            seq_records = future.result()
            for record in seq_records:
                yield record

#########################################
# Sauvegarde en batch via multiprocessing
#########################################

def save_mutated_transcripts_in_batches_parallel(transcripts_fasta, mutations_df, GTF_exons, mutated_output_fasta, batch_size=100, max_workers=32):
    batch = []
    with open(mutated_output_fasta, "w") as mutated_out:
        for mutated_record in _generate_mutated_transcripts_parallel(transcripts_fasta, mutations_df, GTF_exons, max_workers=max_workers):
            batch.append(mutated_record)
            if len(batch) >= batch_size:
                SeqIO.write(batch, mutated_out, "fasta")
                batch = []
        if batch:
            SeqIO.write(batch, mutated_out, "fasta")

####################################
# Sauvegarde des transcrits combinés
####################################

def save_combined_transcripts(transcripts_fasta, mutated_output_fasta, combined_output_fasta):
    original_transcripts = list(SeqIO.parse(transcripts_fasta, "fasta"))
    mutated_transcripts = list(SeqIO.parse(mutated_output_fasta, "fasta"))
    combined_transcripts = original_transcripts + mutated_transcripts
    with open(combined_output_fasta, "w") as combined_out:
        SeqIO.write(combined_transcripts, combined_out, "fasta")

####################
# Pipeline principal
####################

def main():
    exon_parquet = snakemake.input.exon_parquet 
    transcripts_fasta = snakemake.input.transcripts_fasta
    vcf_file = snakemake.input.vcf
    mutated_output_fasta = snakemake.output.mutated_output_fasta
    combined_output_fasta = snakemake.output.combined_output_fasta  
    
    try:
        vcf_df = read_vcf(vcf_file)
        GTF_exons = load_exons(exon_parquet)
        save_mutated_transcripts_in_batches_parallel(transcripts_fasta, vcf_df, GTF_exons, mutated_output_fasta, batch_size=100, max_workers=32)
        save_combined_transcripts(transcripts_fasta, mutated_output_fasta, combined_output_fasta)
    except Exception as e:
        logging.critical(f"Erreur critique dans le pipeline principal: {e}")

if __name__ == "__main__":
    
    main()

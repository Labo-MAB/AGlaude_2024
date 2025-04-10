import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import logging
from itertools import chain, combinations
import concurrent.futures

# Configuration du logger
logging.basicConfig(
    filename="logs/apply_variants.log",
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Variables globales pour multiprocessing
global_genome = None
global_gtf_exons = None

def init_worker(fasta_path, exon_parquet_path):
    global global_genome, global_gtf_exons
    global_genome = {
        record.id.split('.')[0]: record
        for record in SeqIO.parse(fasta_path, "fasta")
    }
    global_gtf_exons = pd.read_parquet(exon_parquet_path)
    global_gtf_exons['exon_interval'] = pd.IntervalIndex.from_arrays(
        global_gtf_exons['start'], global_gtf_exons['end'], closed="both"
    )

def read_vcf(vcf_file):
    mutations = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('##') or line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            mutations.append((fields[0], int(fields[1]), fields[3], fields[4]))
    mutations_df = pd.DataFrame(mutations, columns=['chromosome', 'position', 'ref_nucleotide', 'alt_nucleotide'])
    mutations_df['position'] = mutations_df['position'].astype(np.uint32)
    logging.info("Fichier VCF chargé avec succès.")
    return mutations_df

def all_subsets(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

def _print_transcript_exons(transcript_id):
    global global_gtf_exons, global_genome
    exons_info = global_gtf_exons[global_gtf_exons['transcript_id'] == transcript_id]
    if exons_info.empty or transcript_id not in global_genome:
        logging.warning(f"Transcript {transcript_id} introuvable ou sans exon.")
        return None, None

    strand = exons_info.iloc[0]['strand']
    sorted_exons = exons_info.sort_values(by='start', ascending=(strand == '+'))
    transcript_seq = str(global_genome[transcript_id].seq)
    genomic_to_transcript = {}
    reconstructed = []
    cumulative = 0

    for _, exon in sorted_exons.iterrows():
        exon_len = exon['end'] - exon['start'] + 1
        positions = range(exon['start'], exon['end'] + 1) if strand == '+' else range(exon['end'], exon['start'] - 1, -1)
        for j, genomic_pos in enumerate(positions):
            genomic_to_transcript[genomic_pos] = cumulative + j
        exon_seq = transcript_seq[cumulative: cumulative + exon_len]
        reconstructed.append(exon_seq)
        cumulative += exon_len

    reconstructed_seq = "".join(reconstructed)
    return reconstructed_seq, genomic_to_transcript

def _apply_mutations_to_transcript(transcript_id, transcript_seq, genomic_to_transcript, mutations):
    mutation_list = [m for m in mutations if m['position'] in genomic_to_transcript]
    for m in mutation_list:
        m['transcript_pos'] = genomic_to_transcript[m['position']]
    mutation_list.sort(key=lambda x: x['transcript_pos'])

    offset = 0
    new_seq = transcript_seq
    for m in mutation_list:
        pos = m['transcript_pos'] + offset
        ref, alt = m['ref_nucleotide'], m['alt_nucleotide']
        if new_seq[pos: pos + len(ref)] != ref:
            if len(ref) <= 5:
                logging.warning(
                    f"MUTATION_FAIL | {transcript_id} | pos={pos} | ref={ref} | found={new_seq[pos: pos + len(ref)]} | alt={alt}"
                )
            continue
        new_seq = new_seq[:pos] + alt + new_seq[pos+len(ref):]
        offset += len(alt) - len(ref)
    return new_seq

def _process_transcript_worker(transcript_id, mutations):
    global global_gtf_exons, global_genome
    logging.info(f"Début du traitement pour {transcript_id} ({len(mutations)} mutations)")
    exons_info = global_gtf_exons[global_gtf_exons['transcript_id'] == transcript_id]
    if exons_info.empty or exons_info.iloc[0]['gene_biotype'] == 'lncRNA':
        return []

    transcript_seq, genomic_to_transcript = _print_transcript_exons(transcript_id)
    if transcript_seq is None:
        logging.error(f"Impossible de reconstruire le transcript {transcript_id}.")
        return []

    seq_records = []
    for comb in all_subsets(mutations):
        comb_list = list(comb)
        if any(m['position'] not in genomic_to_transcript for m in comb_list):
            continue
        new_seq = _apply_mutations_to_transcript(transcript_id, transcript_seq, genomic_to_transcript, comb_list)
        signature = "_".join([f"{m['ref_nucleotide']}({m['alt_nucleotide']})" for m in comb_list])
        new_id = f"{transcript_id}_mut_({signature})"
        new_description = global_genome[transcript_id].description
        seq_record = SeqRecord(Seq(new_seq), id=new_id, description=new_description)
        seq_records.append(seq_record)
    return seq_records

def _generate_mutated_transcripts_parallel(transcripts_fasta, mutations_df, exon_parquet, max_workers=32):
    GTF_exons = pd.read_parquet(exon_parquet)
    transcript_mutations = {}
    for _, mutation in mutations_df.iterrows():
        chromosome, position = mutation['chromosome'], mutation['position']
        matching_exons = GTF_exons[(GTF_exons['chromosome'] == chromosome) & (GTF_exons['start'] <= position) & (GTF_exons['end'] >= position)]
        for _, exon in matching_exons.iterrows():
            transcript_id = exon['transcript_id']
            transcript_mutations.setdefault(transcript_id, []).append(mutation.to_dict())

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=init_worker,
        initargs=(transcripts_fasta, exon_parquet)
    ) as executor:
        futures = {
            executor.submit(_process_transcript_worker, tid, muts): tid
            for tid, muts in transcript_mutations.items()
        }
        for i, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            logging.info(f"[{i}/{len(futures)}] Transcript traité.")
            seq_records = future.result()
            for record in seq_records:
                yield record

def save_mutated_transcripts_in_batches_parallel(transcripts_fasta, mutations_df, exon_parquet, mutated_output_fasta, batch_size=100, max_workers=32):
    batch = []
    with open(mutated_output_fasta, "w") as mutated_out:
        for mutated_record in _generate_mutated_transcripts_parallel(transcripts_fasta, mutations_df, exon_parquet, max_workers=max_workers):
            batch.append(mutated_record)
            if len(batch) >= batch_size:
                SeqIO.write(batch, mutated_out, "fasta")
                batch = []
        if batch:
            SeqIO.write(batch, mutated_out, "fasta")

def main():
    exon_parquet = snakemake.input.exon_parquet
    transcripts_fasta = snakemake.input.transcripts_fasta
    vcf_file = snakemake.input.vcf
    mutated_output_fasta = snakemake.output.mutated_output_fasta
    try:
        vcf_df = read_vcf(vcf_file)
        save_mutated_transcripts_in_batches_parallel(
            transcripts_fasta, vcf_df, exon_parquet,
            mutated_output_fasta, batch_size=100, max_workers=32
        )
        logging.info("Pipeline terminé avec succès.")
    except Exception as e:
        logging.critical(f"Erreur critique dans le pipeline principal: {e}")

if __name__ == "__main__":
    main()

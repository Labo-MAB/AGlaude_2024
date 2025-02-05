import os  
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from colorama import init, Fore, Style # a supprimer 


def read_vcf(vcf_file):
    mutations = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('##') or line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chromosome = fields[0]
            position = int(fields[1])
            ref_nucleotide = fields[3]
            alt_nucleotide = fields[4]
            info_field = fields[7]
            gene_name = '.' 

            for item in info_field.split(';'):
                if item.startswith('ANN='):
                    gene_info = item.split('|')
                    if len(gene_info) > 3:
                        gene_name = gene_info[3]
            mutations.append((chromosome, position, ref_nucleotide, alt_nucleotide, gene_name))
    
    mutations_df = pd.DataFrame(mutations, columns=['chromosome', 'position', 'ref_nucleotide', 'alt_nucleotide', 'gene_name'])
    print("Mutations DataFrame:\n", mutations_df) 
    return mutations_df

def find_exons_for_mutations(vcf_df, gtf_file):
    exons = []
    for _, mutation in vcf_df.iterrows():
        found_exon = False
        with open(gtf_file, 'r') as gtf:
            for line in gtf:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if fields[2] == 'exon':
                    chromosome = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    attributes = fields[8].split(';')
                    gene_name = None
                    transcript_id = None  # Ajout de transcript_id

                    for attribute in attributes:
                        if 'gene_name' in attribute:
                            gene_name = attribute.split(' ')[-1].replace('"', '')
                        elif 'transcript_id' in attribute:  # Récupérer le transcript_id
                            transcript_id = attribute.split(' ')[-1].replace('"', '')

                    if (mutation['chromosome'] == chromosome and 
                        start <= mutation['position'] <= end and 
                        mutation['gene_name'] == gene_name):
                        exons.append({
                            'chromosome': chromosome,
                            'position': mutation['position'],
                            'ref_nucleotide': mutation['ref_nucleotide'],
                            'alt_nucleotide': mutation['alt_nucleotide'],
                            'gene_name': gene_name,
                            'transcript_id': transcript_id,  # Ajouter le transcript_id
                            'exon_start': start,
                            'exon_end': end,
                            'strand': strand
                        })
                        found_exon = True
                        break
        
        if not found_exon:
            print(f"Warning: Mutation at position {mutation['position']} not found in any exon for gene {mutation['gene_name']}.")
    
    exons_df = pd.DataFrame(exons)
    print("Annotated Exons DataFrame:\n", exons_df)
    return exons_df


def extract_genomic_sequence(genome_file, exons_df):
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    mutated_sequences = []

    #colorama
    init(autoreset=True)

    for _, exon in exons_df.iterrows():
        chromosome = exon['chromosome']
        position = exon['position']
        alt_nucleotide = exon['alt_nucleotide']
        strand = exon['strand']
        gene_name = exon['gene_name']
        transcript_id = exon['transcript_id'] 

        chromosome_key = f"chr{chromosome}" if f"chr{chromosome}" in genome else chromosome
        
        if chromosome_key in genome:
            seq_record = genome[chromosome_key]
            start_index = exon['exon_start'] - 1  # 0-indexed
            end_index = exon['exon_end']  #  fin exclusive

            if start_index < 0 or end_index > len(seq_record.seq):
                print(f"{Fore.YELLOW}Warning: Indices out of range for chromosome {chromosome} at position {position}.{Style.RESET_ALL}")
                continue

            genomic_sequence = str(seq_record.seq[start_index:end_index])
            if strand == '-':
                genomic_sequence = str(seq_record.seq[start_index:end_index].reverse_complement())

            ref_sequence_start = exon['position'] - 1 - start_index
            ref_sequence_end = ref_sequence_start + len(exon['ref_nucleotide'])

            if 0 <= ref_sequence_start < len(genomic_sequence) and ref_sequence_end <= len(genomic_sequence):
                ref_sequence = genomic_sequence[ref_sequence_start:ref_sequence_end]
                if ref_sequence == exon['ref_nucleotide']:
                    modified_sequence = (
                        genomic_sequence[:ref_sequence_start] +
                        alt_nucleotide +
                        genomic_sequence[ref_sequence_end:]
                    )
                    mutated_sequences.append({
                        'chromosome': chromosome,
                        'position': position,
                        'strand': strand,
                        'original_sequence': genomic_sequence,
                        'modified_sequence': modified_sequence,
                        'gene_name': gene_name,
                        'transcript_id': transcript_id
                    })
                else:
                    print(f"Extracted sequence: {genomic_sequence}")
                    print(f"{Fore.YELLOW}Warning: Reference nucleotide mismatch at {position}. Found: {ref_sequence}, Expected: {exon['ref_nucleotide']}{Style.RESET_ALL}")
            else:
                print(f"{Fore.YELLOW}Warning: Indices out of range for chromosome {chromosome} at position {position}.{Style.RESET_ALL}")
        else:
            print(f"{Fore.YELLOW}Warning: Chromosome {chromosome} not found in the genome.{Style.RESET_ALL}")

    for mutated in mutated_sequences:
        print(f"Chromosome: {mutated['chromosome']}, Position: {mutated['position']}, Strand: {mutated['strand']}")

        original_sequence = mutated['original_sequence']
        modified_sequence = mutated['modified_sequence']

        highlighted_modified_sequence = ''.join(
            f"{Fore.RED}{char}{Style.RESET_ALL}" 
            if i >= len(original_sequence) or char != original_sequence[i] else char
            for i, char in enumerate(modified_sequence)
        )

        print(f"Original Sequence: {original_sequence}")
        print(f"Modified Sequence: {highlighted_modified_sequence}")
        print(f"Gene Name: {mutated['gene_name']}\n")
    
    mutated_sequences_df = pd.DataFrame(mutated_sequences)
    return mutated_sequences_df




def save_mutated_transcripts(transcripts_fasta, mutated_sequences_df, mutated_output_fasta, combined_output_fasta):

    transcripts = {record.id: record for record in SeqIO.parse(transcripts_fasta, "fasta")}
    combined_transcripts = list(transcripts.values())  
    mutated_transcripts = []  

    grouped_mutations = mutated_sequences_df.groupby("transcript_id")

    for transcript_id, group in grouped_mutations:
        if transcript_id not in transcripts:
            print(f"Warning: Transcript ID {transcript_id} not found in transcriptome.")
            continue

        transcript_record = transcripts[transcript_id]
        transcript_seq = str(transcript_record.seq)

        # Appliquer les mutations successivement
        updated_transcript_seq = transcript_seq
        mutation_descriptions = []

        for _, mutation in group.iterrows():
            original_sequence = mutation['original_sequence']
            modified_sequence = mutation['modified_sequence']

            start_pos = updated_transcript_seq.find(original_sequence)
            if start_pos == -1:
                print(f"Warning: Original sequence not found in transcript for {transcript_id}.")
                continue

            end_pos = start_pos + len(original_sequence)

            diff = [
                (i, orig, mod) for i, (orig, mod) in enumerate(zip(original_sequence, modified_sequence))
                if orig != mod
            ]

            for diff_index, original_nuc, modified_nuc in diff:
                mutation_pos_in_transcript = start_pos + diff_index
                mutation_descriptions.append(f"pos {mutation_pos_in_transcript}: {original_nuc}>{modified_nuc}")

                updated_transcript_seq = (
                    updated_transcript_seq[:mutation_pos_in_transcript] +
                    modified_nuc +
                    updated_transcript_seq[mutation_pos_in_transcript + 1:]
                )
        mutation_summary = "; ".join(mutation_descriptions)
        mutated_transcript_record = SeqRecord(
            Seq(updated_transcript_seq),
            id=f"{transcript_id}_mutated",
            description=f"{transcript_record.description} | mutations: {mutation_summary}"
        )

        mutated_transcripts.append(mutated_transcript_record)
        combined_transcripts.append(mutated_transcript_record)

    with open(mutated_output_fasta, "w") as mutated_out:
        SeqIO.write(mutated_transcripts, mutated_out, "fasta")

    with open(combined_output_fasta, "w") as combined_out:
        SeqIO.write(combined_transcripts, combined_out, "fasta")
        
    print(f"Mutated transcripts saved to: {mutated_output_fasta}")
    print(f"Combined transcriptome saved to: {combined_output_fasta}")
    
def main():
    #fasta_file = snakemake.input.fasta
    #vcf_file = snakemake.input.vcf
    #gtf_file = snakemake.input.gtf
    #output_transcriptome_file = snakemake.input.transcriptome
    #output_file = snakemake.output.fasta
    fasta_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/data/references/genome_fa/homo_sapiens_genome.fa"
    vcf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/results/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant_tronque.vcf"
    gtf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/data/references/gtf/homo_sapiens.gtf"
    output_file = "/mnt/c/Users/Antho/Documents/nouveau_projet/workflow/results/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/transcrits_variants.fa"
    mutated_output_fasta = "/mnt/c/Users/Antho/Documents/nouveau_projet/workflow/results/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/mutated_transcripts.fa"
    combined_output_fasta = "/mnt/c/Users/Antho/Documents/nouveau_projet/workflow/results/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/combined_transcripts.fa"
    transcripts_fasta = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/data/references/transcriptome.fa"
    vcf_df = read_vcf(vcf_file)
    annotated_exons_df = find_exons_for_mutations(vcf_df, gtf_file)
    mutated_sequences_df = extract_genomic_sequence(fasta_file, annotated_exons_df)
    save_mutated_transcripts(transcripts_fasta, mutated_sequences_df, mutated_output_fasta, combined_output_fasta)
if __name__ == "__main__":
    main()

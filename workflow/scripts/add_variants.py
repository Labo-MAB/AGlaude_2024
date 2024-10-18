import sys
from Bio import SeqIO
import vcf

# Variables passées depuis Snakemake
vcf_file = sys.argv[1]  # Le chemin vers le fichier VCF
abundance_tsv = sys.argv[2]  # Le chemin vers le fichier TSV (abondances)
output_fasta = sys.argv[3]  # Le chemin vers la sortie FASTA

# Charger les transcrits de référence depuis un fichier FASTA (doit être fourni)
reference_fasta = "path/to/transcripts_reference.fasta"
reference_transcripts = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))

# Charger le fichier VCF avec les variants
vcf_reader = vcf.Reader(open(vcf_file, 'r'))

# Appliquer les variants sur les séquences
for record in vcf_reader:
    chrom = record.CHROM  # Chromosome ou identifiant du transcrit
    pos = record.POS      # Position du variant
    ref = record.REF      # Allèle de référence
    alt = record.ALT      # Allèle alternatif (mutation)
    
    if chrom in reference_transcripts:
        seq = reference_transcripts[chrom].seq
        mutated_seq = seq[:pos-1] + alt[0] + seq[pos:]  # Appliquer la mutation
        reference_transcripts[chrom].seq = mutated_seq  # Mettre à jour la séquence

# Enregistrer le fichier FASTA avec les variants
with open(output_fasta, "w") as fasta_out:
    SeqIO.write(reference_transcripts.values(), fasta_out, "fasta")

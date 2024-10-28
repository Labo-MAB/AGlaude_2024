import sys
from Bio import SeqIO

# Accéder aux fichiers d'entrée et de sortie via Snakemake
fasta_file = snakemake.input.fasta  # Fichier FASTA des transcrits
vcf_file = snakemake.input.vcf       # Fichier VCF filtré
gtf_file = snakemake.input.gtf       # Fichier GTF

# Fichier de sortie
output_file = snakemake.output[0]    # Fichier de sortie pour les transcrits avec mutations

# Logique pour appliquer les mutations aux transcrits
# Exemple simple de lecture du fichier FASTA
def read_fasta(file):
    with open(file, "r") as handle:
        return {record.id: str(record.seq) for record in SeqIO.parse(handle, "fasta")}

# Exemple simple de lecture du fichier VCF pour obtenir les mutations
def read_vcf(file):
    mutations = []
    with open(file, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue  # Ignorer les lignes de header
            fields = line.strip().split("\t")
            # Supposons que la mutation est dans la colonne 4 (allele d'alternative)
            mutations.append((fields[0], int(fields[1]), fields[4]))  # (chromosome, position, mutation)
    return mutations

# Appliquer les mutations aux séquences
def apply_mutations(transcripts, mutations):
    for chromosome, position, mutation in mutations:
        # Appliquer la mutation au bon transcrit ici
        # Cela nécessite de mapper la position au transcrit
        pass  # Logique d'application de mutation à implémenter

# Lire les fichiers
transcripts = read_fasta(fasta_file)
mutations = read_vcf(vcf_file)

# Appliquer les mutations
apply_mutations(transcripts, mutations)

# Écrire le fichier de sortie avec les transcrits d'origine et les versions mutées
with open(output_file, "w") as out_handle:
    for transcript_id, seq in transcripts.items():
        out_handle.write(f">{transcript_id}\n{seq}\n")

import sys
import pandas as pd
from Bio import SeqIO

fasta_in = sys.argv[1]
gtf_parquet = sys.argv[2]
output_dir = sys.argv[3]

# Charger l'exon dataframe
exon_df = pd.read_parquet(gtf_parquet)

# Créer un mapping transcript -> chrom
transcript_chrom = exon_df.drop_duplicates("transcript_id")[["transcript_id", "chromosome"]]
transcript_chrom = transcript_chrom.rename(columns={"chromosome": "chrom"})  # flemme

# Calculer dynamiquement tous les chromosomes pour "rest"
all_chroms = set(transcript_chrom["chrom"].unique())
excluded = set([str(c) for c in range(1, 19)])
rest_chroms = list(all_chroms - excluded)
rest_chroms.sort()  

# Tu peux afficher pour vérifier la première fois
print(f"Chromosomes dans rest: {rest_chroms}")

GROUPS = {
    "1_to_2": ["1", "2"],
    "3_to_4": ["3", "4"],
    "5_to_7": ["5", "6", "7"],
    "7_to_10": ["7", "8", "9", "10"],
    "10_to_13": ["10", "11", "12", "13"],
    "14_to_18": ["14", "15", "16", "17", "18"],
    "rest": rest_chroms  
}

# Construire un dict groupe -> set(transcripts)
group_to_transcripts = {}
for group, chroms in GROUPS.items():
    subset = transcript_chrom[transcript_chrom["chrom"].isin(chroms)]
    group_to_transcripts[group] = set(subset["transcript_id"])

# Lire le FASTA et écrire pour chaque groupe
records_by_group = {group: [] for group in GROUPS}

for record in SeqIO.parse(fasta_in, "fasta"):
    transcript_id = record.id.split()[0]
    for group, transcripts in group_to_transcripts.items():
        if transcript_id in transcripts:
            records_by_group[group].append(record)

for group, records in records_by_group.items():
    output_path = f"{output_dir}/{group}_transcripts.fa"
    with open(output_path, "w") as out_f:
        SeqIO.write(records, out_f, "fasta")

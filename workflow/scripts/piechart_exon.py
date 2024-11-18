import pandas as pd
import matplotlib.pyplot as plt

import io

# Charger le fichier VCF (ignorer les lignes de commentaires commençant par ##)
vcf_file = "F:/breast_cancer/workflow/results/fraction/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant.vcf"
with open(vcf_file, 'r') as f:
    lines = [line for line in f if not line.startswith('##')]

df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')

# Renommer les colonnes pour plus de clarté
df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

# Extraire les informations des annotations (champ INFO, clé ANN)
annotations = []
for info in df["INFO"]:
    ann_field = [field for field in info.split(';') if field.startswith('ANN=')]
    if ann_field:
        ann_value = ann_field[0].split('=')[1]
        annotations.append(ann_value)
    else:
        annotations.append("")

# Extraire les types d'annotations
mutation_types = []
for ann in annotations:
    if ann:
        mutation = ann.split('|')[1]  # Le type de région est dans la 2e colonne après '|'
        mutation_types.append(mutation)
    else:
        mutation_types.append("Unknown")

# Ajouter les types d'annotations au DataFrame
df["Mutation_Type"] = mutation_types

# Compter les occurrences des types de mutation
mutation_counts = df["Mutation_Type"].value_counts()

# Créer un graphique en camembert
plt.figure(figsize=(8, 8))
mutation_counts.plot.pie(autopct='%1.1f%%', startangle=90, colormap='tab10')
plt.title("Répartition des mutations par type de région")
plt.ylabel("")  # Supprimer l'étiquette de l'axe Y
print("bob")
plt.show()

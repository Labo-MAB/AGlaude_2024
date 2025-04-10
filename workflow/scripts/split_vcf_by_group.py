#!/usr/bin/env python3

import os
import subprocess

# Entrées et sorties fournies automatiquement par Snakemake
vcf_file = snakemake.input.vcf
output_path = snakemake.output.vcf
chromosomes = snakemake.params.chromosomes
include_rest = snakemake.params.get("include_rest", False)

# S'assurer que le dossier de sortie existe
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# Extraire les chromosomes présents dans le VCF
all_chroms = subprocess.check_output(
    f"bcftools query -f '%CHROM\\n' {vcf_file} | sort -u", shell=True
).decode().split()

# Ajout automatique des chromosomes non listés
explicit_chroms = set(chromosomes)
if include_rest:
    special_chroms = [c for c in all_chroms if c not in explicit_chroms]
    chromosomes += special_chroms

# Extraire l'en-tête du VCF
header = subprocess.check_output(f"grep '^#' {vcf_file}", shell=True).decode()

# Préparer la commande pour extraire les chromosomes
chrom_args = ' '.join([f"-r {c}" for c in chromosomes])
cmd = f"bcftools view {chrom_args} {vcf_file} | grep -v '^#'"
body = subprocess.check_output(cmd, shell=True).decode()

# Écrire le fichier final
with open(output_path, "w") as f:
    f.write(header)
    f.write(body)

print(f"[Split] {output_path} écrit avec {len(chromosomes)} chromosomes.")

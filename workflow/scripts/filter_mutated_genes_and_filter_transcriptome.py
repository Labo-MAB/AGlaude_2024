import pandas as pd
import sys

import re

def extract_mutated_genes_and_filter_transcriptome(vcf_file, transcriptome_file, output_file):
    mutations = []

    # Lire le fichier VCF et extraire les informations des mutations
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                continue  # Ignorer les lignes de header
            if line.startswith('#'):
                continue  # Ignorer la ligne d'en-tête des colonnes

            fields = line.strip().split('\t')
            info_field = fields[7]
            gene_name = '.'
            is_intergenic = False

            # Vérification des annotations des gènes
            for item in info_field.split(';'):
                if item.startswith('ANN='):
                    gene_info = item.split('|')
                    if len(gene_info) > 3:
                        gene_name = gene_info[3]  # Extraire le nom du gène
                    if 'intergenic_region' in gene_info[1]:
                        is_intergenic = True
                    break  # On s'arrête dès que l'annotation du gène est trouvée

            # Si ce n'est pas un gène intergénique, on ajoute le gène à la liste des mutations
            if not is_intergenic and gene_name != '.':
                mutations.append(gene_name)

    mutated_genes = set(mutations)
    print("Gènes mutés : ", mutated_genes)

    # Filtrer le transcriptome en utilisant la liste des gènes mutés
    with open(transcriptome_file, 'r') as transcriptome, open(output_file, 'w') as output:
        write_transcript = False
        for line in transcriptome:
            if line.startswith(">"):  # Début d'un transcrit
                transcript_name = line.strip().lstrip(">")
                
                # Extraire le nom du gène depuis l'en-tête (par exemple, 'gene=PI4KB')
                match = re.search(r"gene=([^ ]+)", line)
                if match:
                    gene_name_in_header = match.group(1)
                    if gene_name_in_header in mutated_genes:
                        write_transcript = True
                    else:
                        write_transcript = False
                else:
                    write_transcript = False  # Si aucun nom de gène trouvé, ne pas écrire

            if write_transcript:
                output.write(line)

def main():
    transcriptome_file = snakemake.input.transcriptome
    vcf_file = snakemake.input.vcf
    output_file = snakemake.output.filtered_transcriptome
    extract_mutated_genes_and_filter_transcriptome(vcf_file, transcriptome_file, output_file)

if __name__ == "__main__":
    main()

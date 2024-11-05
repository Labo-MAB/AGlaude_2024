import pandas as pd
from Bio import SeqIO

def read_vcf(vcf_file):
    mutations = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                continue  # Ignorer les lignes de commentaires
            if line.startswith('#'):
                continue  # Ignorer l'en-tête
            
            fields = line.strip().split('\t')
            
            # Extraire les informations nécessaires
            chromosome = fields[0]
            position = int(fields[1])
            ref_nucleotide = fields[3]
            alt_nucleotide = fields[4]
            info_field = fields[7]
            
            # Extraire le nom du gène ou la région depuis INFO
            gene_name = '.'
            is_intergenic = False
            for item in info_field.split(';'):
                if item.startswith('ANN='):
                    # Extraire le nom du gène après "ANN=" et le premier '|' (il peut y avoir d'autres annotations)
                    gene_info = item.split('|')
                    if len(gene_info) > 3:
                        gene_name = gene_info[3]
                    # Vérifier si la région est intergénique
                    if 'intergenic_region' in gene_info:
                        is_intergenic = True
                    break
            
            # Ajouter l'entrée à la liste des mutations
            mutations.append((chromosome, position, ref_nucleotide, alt_nucleotide, gene_name, is_intergenic))
    
    # Convertir la liste en DataFrame
    mutations_df = pd.DataFrame(mutations, columns=['chromosome', 'position', 'ref_nucleotide', 'alt_nucleotide', 'gene_name', 'is_intergenic'])
    pd.set_option('display.max_rows', None)  # Affiche toutes les lignes
    pd.set_option('display.max_columns', None)  # Affiche toutes les colonnes

    return mutations_df

def read_gtf(gtf_file):
    gene_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue  # Ignorer les lignes de commentaire
            fields = line.strip().split('\t')
            if fields[2] == 'gene':  # On s'intéresse aux gènes
                # Extraire le nom du gène (gene_name) et non l'identifiant génomique
                attributes = fields[8].split(';')
                gene_name = None
                for attribute in attributes:
                    if 'gene_name' in attribute:
                        gene_name = attribute.split(' ')[-1].replace('"', '')
                        break
                if gene_name:
                    gene_dict[gene_name] = {
                        'chromosome': fields[0],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6]
                    }
    gene_name = "LINC01661"
    if gene_name in gene_dict:
        print(f"{gene_name} trouvé dans gene_dict :")
        print(gene_dict[gene_name])
    else:
        print(f"{gene_name} n'est pas trouvé dans gene_dict.")
#    print(gene_dict)
    return gene_dict


def get_sequence_from_genome(genome_file, chromosome, start, end):
    sequence = ''
    for record in SeqIO.parse(genome_file, "fasta"):
        if record.id == chromosome:
            sequence = str(record.seq[start-1:end])  # Ajuster pour la position 0
            break
    return sequence
def modify_transcript(vcf_file, gtf_file, fasta_file, output_file):
    mutations_df = read_vcf(vcf_file)
    gene_dict = read_gtf(gtf_file)

    with open(output_file, 'w') as out:
        out.write("Le fichier a bien été généré avec les mutations suivantes :\n\n")

        for index, row in mutations_df.iterrows():
            # Vérifier si la mutation est intergénique
            if row.get('is_intergenic', False):  # Par défaut, considère comme False si la colonne n'existe pas
                continue  # Passer à la mutation suivante si is_intergenic est True

            chromosome = row['chromosome']
            position = row['position']
            ref_nucleotide = row['ref_nucleotide']
            alt_nucleotide = row['alt_nucleotide']
            gene_name = row['gene_name']

            if gene_name in gene_dict:
                gene_info = gene_dict[gene_name]
                start = gene_info['start']
                end = gene_info['end']
                
                # Vérifier si la position est dans les limites du gène
                if start <= position <= end:
                    sequence = get_sequence_from_genome(fasta_file, chromosome, start, end)
                    
                    # Calculer l'index et vérifier qu'il est dans les limites de la séquence
                    index_in_sequence = position - start
                    if 0 <= index_in_sequence < len(sequence):
                        # Vérification du nucléotide de référence
                        if sequence[index_in_sequence] == ref_nucleotide:
                            # Appliquer la mutation
                            mutated_sequence = sequence[:index_in_sequence] + alt_nucleotide + sequence[index_in_sequence + 1:]
                            # Écrire dans le fichier de sortie
                            out.write(f">{gene_name}_pos{position}\n{mutated_sequence}\n")
                        else:
                            print(f"Warning: Reference nucleotide does not match at {chromosome}:{position}")
                    else:
                        print(f"Warning: Calculated index {index_in_sequence} is out of range for gene {gene_name} on chromosome {chromosome}")
                else:
                    print(f"Warning: Position {position} is out of the range ({start}-{end}) for gene {gene_name}")
            else:
                print(f"Warning: Gene {gene_name} not found in GTF.")

def main():
    fasta_file = snakemake.input.fasta
    vcf_file = snakemake.input.vcf
    gtf_file = snakemake.input.gtf
    output_file = snakemake.output.fasta

    modify_transcript(vcf_file, gtf_file, fasta_file, output_file)

if __name__ == "__main__":
    main()

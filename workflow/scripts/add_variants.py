import pandas as pd
from Bio import SeqIO

def read_vcf(vcf_file):
    mutations = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                continue  
            if line.startswith('#'):
                continue  

            fields = line.strip().split('\t')
            chromosome = fields[0]
            position = int(fields[1])
            ref_nucleotide = fields[3]
            alt_nucleotide = fields[4]
            info_field = fields[7]
            gene_name = '.'
            is_intergenic = False

            for item in info_field.split(';'):
                if item.startswith('ANN='):
                    # Extrait le nom du gène après "ANN=" et le premier '|' 
                    gene_info = item.split('|')
                    if len(gene_info) > 3:
                        gene_name = gene_info[3]
                        print(f"Gene Info 1: {gene_info[1]}, Gene Info 3: {gene_info[3]}")
                    if 'intergenic_region' in gene_info[1]:
                        is_intergenic = True
                    break
            mutations.append((chromosome, position, ref_nucleotide, alt_nucleotide, gene_name, is_intergenic))
    mutations_df = pd.DataFrame(mutations, columns=['chromosome', 'position', 'ref_nucleotide', 'alt_nucleotide', 'gene_name', 'is_intergenic']) # regarder pour classe/array
    return mutations_df

def read_gtf(gtf_file):
    gene_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue  
            fields = line.strip().split('\t')
            if fields[2] == 'gene': 
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
    return gene_dict

def get_sequence_from_genome(genome_file, chromosome, start, end):
    sequence = ''
    chromosome_found = False  # Flag pour savoir si on trouve le chromosome

    # Parcours les records dans le fichier FASTA
    for record in SeqIO.parse(genome_file, "fasta"):
        print(f"Found chromosome: {record.id}")  # Affiche l'ID de chaque chromosome dans le fichier FASTA
        if record.id == str(chromosome):  # Vérifie que l'identifiant du chromosome correspond à celui passé en paramètre
            chromosome_found = True
            # Vérifie que la plage de positions est valide
            if start < 1 or end > len(record.seq):
                print(f"Erreur : La plage {start}-{end} dépasse la longueur de la séquence pour {chromosome} ({len(record.seq)} bases).")
                return sequence  # Retourne une séquence vide si la plage est invalide + message

            sequence = str(record.seq[start-1:end])  # Extraction de la séquence avec ajustement des indices
            break
    if not chromosome_found:
        print(f"Chromosome {chromosome} non trouvé dans le fichier FASTA.")
    elif len(sequence) == 0:
        print(f"Warning: Séquence vide extraite pour {chromosome} de {start} à {end}.")
    return sequence

def modify_transcript(vcf_file, gtf_file, fasta_file, output_file):
    mutations_df = read_vcf(vcf_file)
    gene_dict = read_gtf(gtf_file)

    with open(output_file, 'w') as out:
        out.write("Le fichier a bien été généré avec les mutations suivantes :\n\n")

        for index, row in mutations_df.iterrows():
            if row.get('is_intergenic', False):  # Par défaut, considère comme False si la colonne n'existe pas
                continue  

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
                    print(f"Length of sequence for gene {gene_name} from {start} to {end} on chromosome {chromosome}: {len(sequence)}")
                    
                    # Calculer l'index relatif pour la mutation dans cette séquence
                    index_in_sequence = position - start
                    print(f"Position: {position}, Start: {start}, Index in sequence: {index_in_sequence}")
                    
                    # Vérifier que l'index est dans les limites de la séquence
                    if 0 <= index_in_sequence < len(sequence):
                        # Vérification du nucléotide de référence
                        if sequence[index_in_sequence] == ref_nucleotide:
                            # Appliquer la mutation
                            mutated_sequence = sequence[:index_in_sequence] + alt_nucleotide + sequence[index_in_sequence + 1:]
                            # Écrire dans le fichier de sortie
                            out.write(f">{gene_name}_pos{position}\n{mutated_sequence}\n")
                            print(f"Success: Mutation for {gene_name} at position {position} written to output file.")
                        else:
                            print(f"Warning: Reference nucleotide does not match at {chromosome}:{position}")
                    else:
                        print(f"Warning: Calculated index {index_in_sequence} is out of range for gene {gene_name} at position {position}")
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

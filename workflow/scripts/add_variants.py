import pysam
from Bio import SeqIO
import pandas as my_little_pandas

def modify_transcript(vcf_file, gtf_file, output_fasta):
    with pysam.VariantFile(vcf_file) as vcf, open(gtf_file) as gtf:
        gene_sequences = {}
        for record in SeqIO.parse(gtf, "gtf"):
            gene_name = record.attributes.get("gene_name")
            if gene_name:
                gene_sequences[gene_name] = str(record.seq)

        with open(output_fasta, "w") as fasta_out:
            for variant in vcf:
                pos = variant.pos
                ref = variant.ref
                alt = variant.alts[0] if variant.alts else None
                anno_info = variant.info.get('ANN')

                if anno_info:
                    # Vérifie si une région intergénique
                    annotation_parts = anno_info[0].split('|')
                    if annotation_parts[0] == 'intergenic':
                        print(f"Mutation à la position {pos} ignorée : région intergénique.")
                        continue  # skip si
                    gene_name = annotation_parts[3]

                    if gene_name in gene_sequences:
                        transcript = gene_sequences[gene_name]                       
                        if pos - 1 < len(transcript):  
                            nucleotide = transcript[pos - 1]  # Les positions commencent à 1
                            if nucleotide == ref:
                                mutated_transcript = transcript[:pos - 1] + alt + transcript[pos:]
                                fasta_out.write(f">{pos}_{ref}_{alt}_{gene_name}\n{mutated_transcript}\n")
                            else:
                                print(f"Erreur : Le nucléotide de référence à la position {pos} pour le gène {gene_name} est {nucleotide}, attendu {ref}.")
                        else:
                            print(f"Erreur : La position {pos} est hors limites pour le gène {gene_name}.")
                    else:
                        print(f"Erreur : Le gène {gene_name} n'a pas été trouvé dans le fichier GTF.")
                else:
                    print("Erreur : Pas d'annotation trouvée dans le variant.")

def main():
    fasta_file = snakemake.input.fasta
    vcf_file = snakemake.input.vcf
    gtf_file = snakemake.input.gtf
    output_fasta = snakemake.output.fasta

    modify_transcript(vcf_file, gtf_file, output_fasta)

if __name__ == "__main__":
    main()

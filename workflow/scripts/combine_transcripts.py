import sys
from Bio import SeqIO

def save_combined_transcripts(transcripts_fasta, mutated_output_fastas, combined_output_fasta):
    """
    Crée un transcriptome personnalisé pour un patient en combinant :
    - tous les transcrits originaux (de transcripts_fasta)
    - tous les transcrits mutés générés à partir des différents morceaux de VCF
    """
    original_transcripts = list(SeqIO.parse(transcripts_fasta, "fasta"))
    mutated_transcripts = []
    for mutated_fasta in mutated_output_fastas:
        mutated_transcripts.extend(list(SeqIO.parse(mutated_fasta, "fasta")))
    combined_transcripts = original_transcripts + mutated_transcripts
    with open(combined_output_fasta, "w") as combined_out:
        SeqIO.write(combined_transcripts, combined_out, "fasta")


try:
    # Utilisation directe des noms définis dans la rule
    transcripts_fasta = snakemake.input.transcripts_fasta
    mutated_output_fastas = [
        snakemake.input.mutated_output_fasta1,
        snakemake.input.mutated_output_fasta2,
        snakemake.input.mutated_output_fasta3,
        snakemake.input.mutated_output_fasta4
    ]
    combined_output_fasta = snakemake.output.combined_transcripts

except NameError:
    # Mode test en ligne de commande
    if len(sys.argv) < 3:
        sys.exit("Usage: combine_transcripts.py <transcripts_fasta> <combined_output_fasta> <mutated_output_fastas...>")
    transcripts_fasta = sys.argv[1]
    combined_output_fasta = sys.argv[2]
    mutated_output_fastas = sys.argv[3:]

save_combined_transcripts(transcripts_fasta, mutated_output_fastas, combined_output_fasta)

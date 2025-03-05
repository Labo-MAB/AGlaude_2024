import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import logging

# Configuration du logger
logging.basicConfig(
    filename="logs/apply_variants.log",
    
    level=logging.DEBUG,               
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def read_vcf(vcf_file):
    """Charge le fichier VCF en DataFrame optimisé."""
    mutations = []
    try:
        with open(vcf_file, 'r') as vcf:
            for line in vcf:
                if line.startswith('##') or line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                chromosome = fields[0]
                position = int(fields[1])
                ref_nucleotide = fields[3]
                alt_nucleotide = fields[4]
                info_field = fields[7]
                gene_name = None
                for item in info_field.split(';'):
                    if item.startswith('ANN='):
                        gene_info = item.split('|')
                        if len(gene_info) > 3:
                            gene_name = gene_info[3]
                mutations.append((chromosome, position, ref_nucleotide, alt_nucleotide, gene_name))
        mutations_df = pd.DataFrame(mutations, columns=['chromosome', 'position', 'ref_nucleotide', 'alt_nucleotide', 'gene_name'])
        mutations_df['position'] = mutations_df['position'].astype(np.uint32)
        logging.info("VCF chargé avec succès.")
        print("Mutations DataFrame:\n", mutations_df)
        return mutations_df
    except Exception as e:
        logging.error(f"Erreur lors de la lecture du VCF: {e}")
        raise

def load_exons(exon_parquet):
    """Charge les exons depuis un fichier Parquet pour éviter un recalcul."""
    try:
        GTF_exons = pd.read_parquet(exon_parquet)
        GTF_exons['exon_interval'] = pd.IntervalIndex.from_arrays(GTF_exons['start'], GTF_exons['end'], closed="both")
        logging.info("Exons chargés avec succès.")
        return GTF_exons
    except Exception as e:
        logging.error(f"Erreur lors du chargement des exons: {e}")
        raise

def _print_transcript_exons(transcript_id, GTF_exons, genome):
    """
    Reconstruit la séquence transcriptomique et établit un mapping entre positions génomiques et transcriptomiques.
    
    Retourne :
      - reconstructed_seq : séquence complète du transcript
      - genomic_to_transcript : dictionnaire {position génomique -> position transcriptomique}
    """
    exons_info = GTF_exons[GTF_exons['transcript_id'] == transcript_id]
    if exons_info.empty:
        msg = f"Aucun exon trouvé pour le transcript {transcript_id}."
        logging.warning(msg)
        print(msg)
        return None, None

    strand = exons_info.iloc[0]['strand']
    sorted_exons = exons_info.sort_values(by='start', ascending=(strand == '+'))
    
    try:
        transcript_seq = str(genome[transcript_id].seq)
    except KeyError as e:
        logging.error(f"Transcript {transcript_id} non trouvé dans le FASTA: {e}")
        return None, None

    genomic_to_transcript = {}
    reconstructed = []
    cumulative = 0  # Position dans le transcript reconstruit

    print(f"\n### DEBUG : Reconstruction du transcript {transcript_id}")
    for i, (_, exon) in enumerate(sorted_exons.iterrows(), start=1):
        exon_len = exon['end'] - exon['start'] + 1
        if strand == '+':
            for j, genomic_pos in enumerate(range(exon['start'], exon['end'] + 1)):
                genomic_to_transcript[genomic_pos] = cumulative + j
        else:
            for j, genomic_pos in enumerate(range(exon['end'], exon['start'] - 1, -1)):
                genomic_to_transcript[genomic_pos] = cumulative + j
        exon_seq = transcript_seq[cumulative: cumulative + exon_len]
        reconstructed.append(exon_seq)
        logging.debug(f"Exon {i} - Genomic: {exon['start']}-{exon['end']}, Transcript: {cumulative}-{cumulative+exon_len-1}")
        print(f"Exon {i} - Coordonnées génomiques: {exon['start']} - {exon['end']}, positions transcriptomiques: {cumulative} - {cumulative+exon_len-1}")
        print(f"Séquence exon: {exon_seq}\n")
        cumulative += exon_len

    reconstructed_seq = "".join(reconstructed)
    logging.info(f"Transcript {transcript_id} reconstruit, longueur {len(reconstructed_seq)}.")
    print("Transcript complet reconstruit :", reconstructed_seq)
    print("Longueur =", len(reconstructed_seq))
    return reconstructed_seq, genomic_to_transcript

def apply_mutations_to_transcript(transcript_seq, genomic_to_transcript, mutations):
    """
    Applique une liste de mutations à une séquence transcriptomique.
    """
    mutation_list = []
    for m in mutations:
        genomic_pos = m['position']
        if genomic_pos in genomic_to_transcript:
            m["transcript_pos"] = genomic_to_transcript[genomic_pos]
            mutation_list.append(m)
        else:
            msg = f"Position {genomic_pos} non trouvée dans le mapping pour la mutation {m}."
            logging.warning(msg)
            print(msg)
    
    mutation_list.sort(key=lambda x: x["transcript_pos"])
    offset = 0
    new_seq = transcript_seq
    for m in mutation_list:
        pos = m["transcript_pos"] + offset
        ref = m["ref_nucleotide"]
        alt = m["alt_nucleotide"]
        if new_seq[pos: pos + len(ref)] != ref:
            msg = f"**Mismatch** à la position {pos} dans le transcript: attendu {ref}, trouvé {new_seq[pos: pos + len(ref)]}"
            logging.error(msg)
            print(msg)
        new_seq = new_seq[:pos] + alt + new_seq[pos+len(ref):]
        offset += len(alt) - len(ref)
    return new_seq

def extract_genomic_sequence(transcripts_fasta, mutations_df, GTF_exons):
    """
    Reconstruit la séquence transcriptomique pour chaque transcript concerné et applique toutes les mutations cumulativement.
    Retourne un dictionnaire {transcript_id: {"sequence": séquence mutée, "signature": signature de mutation}}.
    """
    genome = SeqIO.to_dict(SeqIO.parse(transcripts_fasta, "fasta"))
    transcript_mutations = {}  # Regrouper les mutations par transcript

    for idx, mutation in mutations_df.iterrows():
        chromosome = mutation['chromosome']
        position = mutation['position']  # Position génomique (1-indexée)
        ref_nucleotide = mutation['ref_nucleotide']
        alt_nucleotide = mutation['alt_nucleotide']
        gene_name = mutation['gene_name']

        msg = f"Pour la mutation sur le gène {gene_name} à la position {position}: {ref_nucleotide} -> {alt_nucleotide}"
        logging.info(msg)
        print("\n" + msg)

        matching_exons = GTF_exons[
            (GTF_exons['gene_name'] == gene_name) &
            (GTF_exons['chromosome'] == chromosome) &
            (GTF_exons['start'] <= position) &
            (GTF_exons['end'] >= position)
        ]
        logging.debug(f"Matching exons pour {gene_name} à {position}: {matching_exons}")
        print("Matching exons trouvés:")
        print(matching_exons)

        if matching_exons.empty:
            msg = f"Aucune correspondance d'exon pour {gene_name} à la position {position}."
            logging.warning(msg)
            print(msg)
            continue

        for _, exon in matching_exons.iterrows():
            transcript_id = exon['transcript_id']
            if transcript_id not in genome:
                msg = f"Transcript {transcript_id} non trouvé dans le FASTA."
                logging.error(msg)
                print(msg)
                continue
            if transcript_id not in transcript_mutations:
                transcript_mutations[transcript_id] = []
            # On ajoute aussi l'index de la ligne pour identifier précisément la mutation
            mutation_detail = mutation.to_dict()
            mutation_detail['line_index'] = idx
            transcript_mutations[transcript_id].append(mutation_detail)

    mutated_sequences = {}
    for transcript_id, mutations in transcript_mutations.items():
        transcript_seq, genomic_to_transcript = _print_transcript_exons(transcript_id, GTF_exons, genome)
        if transcript_seq is None or genomic_to_transcript is None:
            continue
        try:
            new_seq = apply_mutations_to_transcript(transcript_seq, genomic_to_transcript, mutations)
        except Exception as e:
            logging.error(f"Erreur lors de l'application des mutations pour {transcript_id}: {e}")
            continue
        # Créer une signature de mutation pour renommer le header
        mutations_sorted = sorted(mutations, key=lambda m: m['position'])
        signature_parts = [f"{m['ref_nucleotide']}({m['alt_nucleotide']})" for m in mutations_sorted]
        signature = "_".join(signature_parts)
        mutated_sequences[transcript_id] = {"sequence": new_seq, "signature": signature}
        msg = f"Transcript {transcript_id} : mutations cumulatives appliquées. Longueur originale: {len(transcript_seq)}, longueur modifiée: {len(new_seq)}, signature: {signature}"
        logging.info(msg)
        print("\n" + msg)
    return mutated_sequences

def save_mutated_transcripts(mutated_info, mutated_output_fasta, combined_output_fasta, transcripts_fasta):
    """
    Sauvegarde les transcrits mutés en reprenant le header original et en le modifiant pour inclure le préfixe 'mut_' et la signature.
    Format du nouveau header: >mut_<transcript_id>_(<mutation_signature>)
    """
    original_transcripts = {record.id: record for record in SeqIO.parse(transcripts_fasta, "fasta")}
    mutated_transcripts = []

    for transcript_id, info in mutated_info.items():
        mutated_seq = info["sequence"]
        signature = info["signature"]
        if transcript_id in original_transcripts:
            original_record = original_transcripts[transcript_id]
            new_id = f"mut_{transcript_id}_({signature})"
            new_description = original_record.description  # On garde la description originale
        else:
            new_id = f"mut_{transcript_id}_({signature})"
            new_description = f"Mutated transcript {transcript_id}"
        mutated_record = SeqRecord(Seq(mutated_seq), id=new_id, description=new_description)
        mutated_transcripts.append(mutated_record)

    try:
        with open(mutated_output_fasta, "w") as mutated_out:
            SeqIO.write(mutated_transcripts, mutated_out, "fasta")
        logging.info(f"Transcrits mutés sauvegardés dans {mutated_output_fasta}.")
    except Exception as e:
        logging.error(f"Erreur lors de la sauvegarde des transcrits mutés: {e}")

    try:
        original_list = list(original_transcripts.values())
        combined_transcripts = original_list + mutated_transcripts
        with open(combined_output_fasta, "w") as combined_out:
            SeqIO.write(combined_transcripts, combined_out, "fasta")
        logging.info(f"Transcrits combinés sauvegardés dans {combined_output_fasta}.")
    except Exception as e:
        logging.error(f"Erreur lors de la sauvegarde des transcrits combinés: {e}")

def main():
    exon_parquet = snakemake.input.exon_parquet 
    transcripts_fasta = snakemake.input.transcripts_fasta
    vcf_file = snakemake.input.vcf
    mutated_output_fasta = snakemake.output.mutated_output_fasta
    combined_output_fasta = snakemake.output.combined_output_fasta

    try:
        vcf_df = read_vcf(vcf_file)
        GTF_exons = load_exons(exon_parquet)
        mutated_info = extract_genomic_sequence(transcripts_fasta, vcf_df, GTF_exons)
        save_mutated_transcripts(mutated_info, mutated_output_fasta, combined_output_fasta, transcripts_fasta)
    except Exception as e:
        logging.critical(f"Erreur critique dans le pipeline principal: {e}")

if __name__ == "__main__":
    main()



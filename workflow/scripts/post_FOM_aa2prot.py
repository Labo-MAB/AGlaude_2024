import pickle
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pprint import pprint
from collections import defaultdict
import os

def load_pickle(file_path):
    """Charge le fichier pickle et retourne le dictionnaire associé."""
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    return data

def load_custom_transcriptome(fasta_file):
    """
    Charge le transcriptome personnalisés (contenant WT avec CDS + mutée )
    et retourne un dictionnaire de SeqRecord
    """
    return SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

def group_keys(pickle_data):
    """
    Regroupe les entrées du pickle par transcript_id, séparant 'wt' et 'mut'.
    """
    transcript_entries = {}
    for key, value in pickle_data.items():
        if '_mut_' in key:
            transcript_id = key.split('_mut_')[0]
            if transcript_id not in transcript_entries:
                transcript_entries[transcript_id] = {"wt": [], "mut": []}
            transcript_entries[transcript_id]["mut"].append((key, value))
        else:
            transcript_id = key
            if transcript_id not in transcript_entries:
                transcript_entries[transcript_id] = {"wt": [], "mut": []}
             
#               if transcript_id == "ENST00000832823":

            transcript_entries[transcript_id]["wt"].append((key, value))
            if transcript_id == "ENST00000832823":
                with open("output/debug_ENST00000832823.txt", "a") as debug_f:
                    debug_f.write(f" {transcript_id} | val = {key}\n")

    
    return transcript_entries


def add_ref_cds(valid_dict, transcriptome):
    """
    Pour chaque transcript dans le dictionnaire groupé, ajoute l'information du CDS annoté
    issue du transcriptome (fichier FASTA).

    Exemple : pour un header FASTA tel que 
      >ENST00000641515 gene=OR4F5 CDS=61-1039
    on extrait (61, 1039) et on le stocke dans valid_dict[transcript_id]['ref'].
    """
    
    for transcript_id, data in valid_dict.items():
        if transcript_id in transcriptome:
            record = transcriptome[transcript_id]
            match = re.search(r'CDS=(\d+)-(\d+)', record.description)
            if match:
                start = int(match.group(1))
                end = int(match.group(2))
                data['ref'] = (start, end)
            else:
                data['ref'] = None
        else:
            data['ref'] = None
    return valid_dict

def extract_first_coordinate(entries):
    """
    À partir d'une liste d'entrées (issues du pickle) sous la forme [(clé, valeur), ...],
    retourne le premier tuple trouvé. Si la liste est vide ou mal formée, retourne None.
    """
    if entries:
        val = entries[0][1]
        if val and len(val) > 0:
            return val[0]
    return None

def check_cds_match(wt_coords, ref_coords):
    """
    Vérifie si les coordonnées wt (wild-type, du pickle) correspondent aux coordonnées de référence (du FASTA),
    après avoir converti :
      - wt : base 0 ➔ base 1 (+1)
      - wt : inclut le codon stop ➔ on enlève 3 nt pour comparer avec un FASTA sans stop.
    """
    if not wt_coords and not ref_coords:
        return True
    if not wt_coords or not ref_coords:
        return False

    wt_start, wt_end = wt_coords
    ref_start, ref_end = ref_coords

    # Convertir pickle (FOMOnet)
    wt_start += 1  # 0-based ➔ 1-based
    wt_end += 1
    wt_end -= 3    # enlever le codon stop pour comparer proprement

    start_diff = abs(wt_start - ref_start)
    end_diff = abs(wt_end - ref_end)

    return (start_diff <= 1) and (end_diff <= 2)


def verify_predictions(ref_cds_dict, output_file):
    """
    Vérifie pour chaque transcript si les coordonnées wt issues du pickle correspondent
    aux coordonnées de référence (ref) extraites du FASTA.
    Les transcripts ne passant pas la vérification sont écrits dans output_file et retirés.
    """
    valid_transcripts = {}
    error_entries = []
    for transcript_id, data in ref_cds_dict.items():
        wt_coords = extract_first_coordinate(data.get("wt", []))
        ref_coords = data.get("ref", None)
        if check_cds_match(wt_coords, ref_coords):
            valid_transcripts[transcript_id] = data
        else:
            error_entries.append(f"{transcript_id} | wt: {wt_coords} | ref: {ref_coords}\n")
    with open(output_file, "w") as out:
         out.writelines(error_entries)
    return valid_transcripts





def _parse_mutation_key(mut_key):
    """
    Parse une clé de mutation de forme :
        ENST..._mut_<sig1>_<sig2>_..._pos=<p1>;<p2>;...
    Exemples traités :
        ENST00000473530_mut_G(A)_G(A)_pos=790;1078
        ENST00000699555_mut_C(T)_G(C)_A(G)_A(G)_pos=512;875;984;1056
        ENST00000356006_mut_CGCA(GGCG)_CCCTC(TCCTT)_pos=1181;1205
    Retourne :
        transcript_id (str),
        signature (str) où les mutations sont séparées par ';' (ex. "G(A);T(C)"),
        positions (list[int])
    """
    try:
        # 1) Sépare l’ID de la partie mutations+positions
        transcript_id, rest = mut_key.split('_mut_', 1)
        # 2) Sépare la partie mutations (sig_str) de la partie positions (pos_str)
        sig_str, pos_str = rest.rsplit('_pos=', 1)
    except ValueError:
        raise ValueError(f"Clé mutée non reconnue : {mut_key}")

    # 3) Remplace les underscores par des ';' dans la signature
    #    pour normaliser "G(A)_T(C)_G(A)" → "G(A);T(C);G(A)"
    signature = sig_str.replace('_', ';')

    # 4) Convertit les positions en liste d’entiers
    positions = [int(p) for p in pos_str.split(';')]

    return transcript_id, signature, positions



def _corriger_positions(positions, ref_cds):

    """
    Corrige les positions mutées en les recentrant dans le CDS si elles sont hors bornes.
    Permet de maintenir l’analyse du CDS même en cas de mutations amont ou aval.
    """
    if ref_cds is None:
        return positions  #####   
    ref_start, ref_end = ref_cds
    positions_corrigees = []
    for pos in positions:
        if pos < ref_start:
            positions_corrigees.append(ref_start)
        elif pos > ref_end:
            positions_corrigees.append(ref_end)
        else:
            positions_corrigees.append(pos)
    return positions_corrigees


def analyser_impact(filtered_dict, output_file, transcriptome):
    """
    Analyse les mutations et écrit un fichier TSV listant:
      transcript_id, mutation, positions, corrected_positions, frameshift, in_cds, cds_impact
    Si transcriptome est fourni, compare les codons impactés.
    """
    with open(output_file, "w") as f:
        f.write("transcript_id\tmutation\tpositions\tcorrected_positions\tframeshift\tin_cds\tcds_impact\n")

        for transcript_id, data in filtered_dict.items():
            ref_cds = data.get("ref")
            wt_coords = extract_first_coordinate(data.get("wt", []))

            for mut_key, val in data.get("mut", []):
                try:
                    _, signature, positions = _parse_mutation_key(mut_key)
                    corrected_positions = _corriger_positions(positions, ref_cds)
                    mut_coords = val[0] if val and len(val) > 0 else None

                    # Vérifier frameshift
                    mutation_list = signature.split(";")
                    frameshift = False
                    for ref_alt in mutation_list:
                        ref, alt = ref_alt.split("(")
                        alt = alt.rstrip(")")
                        if len(ref) != len(alt):
                            delta = abs(len(ref) - len(alt))
                            if delta % 3 != 0:
                                frameshift = True
                                break

                    # Vérifier si mutation dans CDS
                    in_cds = any(ref_cds[0] <= pos <= ref_cds[1] for pos in corrected_positions) if ref_cds else False

                    # Déterminer cds_impact
                    if not mut_coords:
                        cds_impact = "cds_lost" if wt_coords else "none"
                    else:
                        if not wt_coords:
                            cds_impact = "cds_gain"
                        else:
                            if not frameshift:
                                if wt_coords == mut_coords:
                                    cds_impact = "none"
                                else:
                                    if not in_cds:
                                        impact_detected = False
                                        for pos in corrected_positions:
                                            if ref_cds and (ref_cds[0] - 6 <= pos <= ref_cds[0] + 4):
                                                cds_impact = "kozak"
                                                impact_detected = True
                                                break
                                        if not impact_detected:
                                            cds_impact = "no_impact"
                                    else:
                                        if transcriptome and transcript_id in transcriptome and ref_cds:
                                            record = transcriptome[transcript_id]
                                            silent = True
                                            mutation_list = signature.split(";")
                                            mut_info = []
                                            for mut_str in mutation_list:
                                                try:
                                                    ref, alt = mut_str.split("(")
                                                    alt = alt.rstrip(")")
                                                    mut_info.append((ref, alt))
                                                except Exception as e:
                                                    print(f"Erreur parsing mutation: {mut_str} ({e})")

                                            for (pos, (ref, alt)) in zip(corrected_positions, mut_info):
                                                rel_pos = pos - ref_cds[0]
                                                codon_start_offset = (rel_pos // 3) * 3
                                                codon_genome_start = ref_cds[0] + codon_start_offset
                                                codon_seq = record.seq[codon_genome_start - 1 : codon_genome_start + 2]
                                                if len(codon_seq) != 3:
                                                    continue
                                                pos_in_codon = pos - codon_genome_start
                                                mut_codon_list = list(str(codon_seq))
                                                mut_codon_list[pos_in_codon] = alt
                                                codon_mut_seq = "".join(mut_codon_list)
                                                aa_wt = Seq(codon_seq).translate()
                                                aa_mut = Seq(codon_mut_seq).translate()
                                                if aa_wt != aa_mut:
                                                    silent = False
                                                    break

                                            cds_impact = "silent" if silent else "missense"
                                        else:
                                            cds_impact = "shift"
                            else:
                                if ref_cds == mut_coords:
                                    cds_impact = "none"
                                else:
                                    wt_len = wt_coords[1] - wt_coords[0] + 1 if wt_coords else 0
                                    mut_len = mut_coords[1] - mut_coords[0] + 1
                                    if mut_len == wt_len:
                                        cds_impact = "shift_with_frameshift"
                                    elif mut_len > wt_len:
                                        cds_impact = "cds_gain"
                                    elif mut_len < wt_len:
                                        cds_impact = "cds_loss"
                                    else:
                                        cds_impact = "unknown"

                    f.write(
                        f"{transcript_id}\t{signature}\t{';'.join(map(str, positions))}\t"
                        f"{';'.join(map(str, corrected_positions))}\t{frameshift}\t{in_cds}\t{cds_impact}\n"
                    )

                except Exception as e:
                    print(f"Erreur analyse {mut_key}: {e}")




#########################################################
# création fasta protein
#########################################################

# --- Fonctions pour la traduction en protéine et l'écriture du FASTA ---


def _parse_mutation_key(mut_key):
    transcript_id, rest = mut_key.split('_mut_', 1)
    sig_str, pos_str = rest.rsplit('_pos=', 1)
    signature = sig_str.replace('_', ';')
    positions = [int(p) for p in pos_str.split(';')]
    return transcript_id, signature, positions

def _get_protein_sequence(record, cds):
    if not cds:
        return None
    start, end = cds
    coding_seq = record.seq[start-1:end]
    protein_seq = coding_seq.translate(to_stop=True)
    return str(protein_seq)

def _get_gene_name(record):
    match = re.search(r'gene=([^\s]+)', record.description)
    return match.group(1) if match else "unknown"

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

def write_fasta_sequences(filtered_dict, transcriptome, output_fasta, diagnostic_tsv):
    proteins_by_transcript = {}
    skipped_infos = []  # ← ici on accumule les lignes du TSV

    skipped_missing = 0
    skipped_no_protein = 0

    for transcript_id, data in filtered_dict.items():
        record = transcriptome.get(transcript_id)
        cds = data.get("ref")

        if record is None:
            print(f"[SKIP] Transcript {transcript_id} absent du transcriptome.")
            skipped_missing += 1
            continue

        # extraire séquence et CDS
        seq_len = len(record.seq)
        cds_len = None
        cause = None

        if not cds:
            cause = "Pas de CDS"
        else:
            start, end = cds
            cds_len = end - start + 1

            if end > seq_len:
                cause = "CDS hors-borne"
            elif cds_len < 3:
                cause = "CDS trop court"
            else:
                cds_seq = record.seq[start-1:end]
                if not set(cds_seq.upper()).issubset({"A", "T", "C", "G"}):
                    cause = "Bases invalides"
                else:
                    try:
                        prot_seq = cds_seq.translate(to_stop=True)
                        if len(prot_seq.strip()) == 0 or str(prot_seq.strip()) == "*":
                            cause = "Traduction vide"
                        else:
                            cause = None
                    except Exception:
                        cause = "Erreur traduction"

        # si une cause a été identifiée → skip + ajouter à diagnostic
        if cause:
            print(f"[SKIP] Pas de protéine WT pour {transcript_id} ({cause})")
            skipped_no_protein += 1
            skipped_infos.append([
                transcript_id,
                f"{cds[0]}-{cds[1]}" if cds else "NA",
                seq_len,
                cds_len if cds_len else "NA",
                cause
            ])
            continue

        # on génère normalement si pas de cause
        prot_seq = cds_seq.translate(to_stop=True)
        prot_seq = str(prot_seq).strip().upper()

        gene_name = _get_gene_name(record)
        proteins_by_transcript[transcript_id] = {}
        proteins_by_transcript[transcript_id][transcript_id] = {
            "sequence": prot_seq,
            "TA": set([transcript_id]),
            "gene_name": gene_name
        }

        for mut_key, val in data.get("mut", []):
            try:
                _, signature, positions = _parse_mutation_key(mut_key)
                seq_mut = list(str(record.seq))
                mutations = signature.split(";")

                for mut_str, pos in zip(mutations, positions):
                    ref_base, alt_base = mut_str.split("(")
                    alt_base = alt_base.rstrip(")")
                    pos_index = pos - 1
                    if pos_index < 0 or pos_index >= len(seq_mut):
                        print(f"[WARN] Pos hors borne {pos} pour {mut_key}")
                        continue
                    seq_mut[pos_index] = alt_base

                if cds:
                    start, end = cds
                    if end > len(seq_mut):
                        print(f"[WARN] CDS end {end} dépasse seq {len(seq_mut)} pour {mut_key}")
                        continue
                    coding_seq = "".join(seq_mut[start - 1:end])
                else:
                    coding_seq = "".join(seq_mut)

                mut_prot_seq = str(Seq(coding_seq).translate(to_stop=True)).strip().upper()

                found = False
                for existing_id, obj in proteins_by_transcript[transcript_id].items():
                    if obj["sequence"] == mut_prot_seq:
                        obj["TA"].add(mut_key)
                        found = True
                        break

                if not found:
                    proteins_by_transcript[transcript_id][mut_key] = {
                        "sequence": mut_prot_seq,
                        "TA": set([mut_key]),
                        "gene_name": gene_name
                    }

            except Exception as e:
                print(f"[ERROR] Mutation {mut_key} : {e}")

    print(f"[WRITE] {len(proteins_by_transcript)} transcripts à écrire")
    print(f"[SKIP] {skipped_missing} transcripts absents du transcriptome")
    print(f"[SKIP] {skipped_no_protein} transcripts sans protéine WT")

    with open(output_fasta, "w") as out_f:
        for transcript_id, subdict in proteins_by_transcript.items():
            for entry_id, obj in subdict.items():
                ta_list = sorted(obj["TA"])
                header = (
                    f">{entry_id}|OS=Homo sapiens|GN={obj['gene_name']}|"
                    f"TA={','.join(ta_list)}|PA="
                )
                out_f.write(header + "\n")
                out_f.write(obj["sequence"] + "\n")

    # Sauvegarde du TSV
    df = pd.DataFrame(skipped_infos, columns=[
        "transcript_id", "CDS", "sequence_length", "CDS_length", "problème"
    ])
    df.to_csv(diagnostic_tsv, sep="\t", index=False)
    print(f"[INFO] Fichier diagnostique écrit dans {diagnostic_tsv}")






def merge_with_uniprot(transcript_fasta, uniprot_fasta, final_output):
    """
    Fusionne les séquences du transcriptome avec UniProt.
    Si une séquence wild-type est identique à une séquence UniProt ou déjà présente localement,
    les mutations silencieuses sont fusionnées dans le champ TA=.
    """

    # Lecture UniProt
    uniprot_records = list(SeqIO.parse(uniprot_fasta, "fasta"))
    uniprot_dict = {str(record.seq).strip(): record for record in uniprot_records}

    # Lecture des transcrits
    transcript_records = list(SeqIO.parse(transcript_fasta, "fasta"))

    # Organisation des séquences par contenu
    transcript_groups = defaultdict(list)
    for record in transcript_records:
        seq_str = str(record.seq).strip()
        transcript_groups[seq_str].append(record)

    final_records = []

    for seq_str, records in transcript_groups.items():
        ids = [rec.id for rec in records]
        mutated_ids = [i for i in ids if "_mut(" in i]
        is_wt_present = any("_mut(" not in i for i in ids)

        if seq_str in uniprot_dict:
            # Correspondance UniProt → fusionne dans TA=
            uniprot_rec = uniprot_dict[seq_str]
            match = re.search(r'(TA=)([^ ]+)', uniprot_rec.description)
            existing_ta = match.group(2).split(',') if match else []

            for mid in mutated_ids:
                if mid not in existing_ta:
                    existing_ta.append(mid)

            new_ta = ",".join(existing_ta)
            new_description = re.sub(r'(TA=)([^ ]+)', f'TA={new_ta}', uniprot_rec.description)
            if 'TA=' not in new_description:
                new_description += f' TA={new_ta}'

            updated_record = SeqRecord(
                Seq(seq_str),
                id=uniprot_rec.id,
                description=new_description
            )
            final_records.append(updated_record)

        elif is_wt_present:
            # Pas UniProt mais déjà un wild-type local → fusionne mutations silencieuses
            ref_record = next((rec for rec in records if "_mut(" not in rec.id), None)
            if ref_record:
                match = re.search(r'(TA=)([^ ]+)', ref_record.description)
                existing_ta = match.group(2).split(',') if match else []

                for mid in mutated_ids:
                    if mid not in existing_ta:
                        existing_ta.append(mid)

                new_ta = ",".join(existing_ta)
                new_description = re.sub(r'(TA=)([^ ]+)', f'TA={new_ta}', ref_record.description)
                if 'TA=' not in new_description:
                    new_description += f' TA={new_ta}'

                ref_record.description = new_description
                final_records.append(ref_record)
        else:
            # Ni UniProt ni wild-type local → ajoute tout tel quel
            for rec in records:
                final_records.append(rec)

    SeqIO.write(final_records, final_output, "fasta")

def diagnostiquer_transcrits_non_traduits(filtered_dict, transcriptome, output_tsv, skipped_ids):
    """
    Produit un fichier TSV listant les transcrits avec CDS mais sans protéine WT traduite,
    avec explication de la cause : CDS hors-borne, trop court, N non traduisibles, etc.
    """
    lignes = []

    for transcript_id in skipped_ids:
        record = transcriptome.get(transcript_id)
        cds = filtered_dict.get(transcript_id, {}).get("ref")

        if not record or not cds:
            continue  # transcript absent ou pas de CDS connu dans filtered_dict

        start, end = cds
        seq_len = len(record.seq)
        cds_len = end - start + 1

        if end > seq_len:
            cause = "CDS hors-borne"
        elif cds_len < 3:
            cause = "CDS trop court"
        else:
            cds_seq = record.seq[start-1:end]
            if not set(cds_seq.upper()).issubset({"A", "T", "C", "G"}):
                cause = "Bases invalides"
            else:
                try:
                    prot = cds_seq.translate(to_stop=True)
                    if len(prot) == 0:
                        cause = "Traduction vide"
                    else:
                        cause = "Inconnue"
                except Exception:
                    cause = "Erreur traduction"

        lignes.append([
            transcript_id,
            f"{start}-{end}",
            seq_len,
            cds_len,
            cause
        ])

    import pandas as pd
    df = pd.DataFrame(
        lignes,
        columns=["transcript_id", "CDS", "sequence_length", "CDS_length", "problème"]
    )
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[INFO] Diagnostic écrit dans : {output_tsv}")



def main():
    pickle_file = "SRR43_tx_orfs.pkl"
    fasta_path = "../FOMOnet/FOMOnet/data/transcriptome.fa"
    fasta_path = "../transcriptome.fa"

#    uniprot_fasta = "/mnt/c/Users/Antho/OneDrive - USherbrooke/Documents/gadph_GTEX/Final_gtex/expression_withcluster/pseudogenes/1_expression_Gtex/references/human-openprot-2_1-refprots+altprots+isoforms-uniprot2022_06_01.fasta"
    uniprot_fasta = "/mnt/c/Users/Anthony/OneDrive - USherbrooke/Documents/gadph_GTEX/Final_gtex/expression_withcluster/pseudogenes/1_expression_Gtex/references/human-openprot-2_1-refprots+altprots+isoforms-uniprot2022_06_01.fasta"

    # Fichiers de sortie
    os.makedirs("SRR43", exist_ok=True)
    final_output_fasta = "SRR43/final_merged.fasta"
    output_file = "SRR43/mauvaise_prediction.txt"
    output_fasta = "SRR43/sequences_proteiques.fasta"
    impact_outfile = "SRR43/mutations_impact.tsv"
    diagnostique_tsv = "SRR43/transcrits_non_traduits.tsv"

    # Chargement des données
    pickle_data = load_pickle(pickle_file)
    transcriptome = load_custom_transcriptome(fasta_path)
    valid_dict = group_keys(pickle_data)                                    # Groupement des transcrits WT / mutés
    ref_cds_dict = add_ref_cds(valid_dict, transcriptome)          # Ajout des coordonnées CDS depuis le fichier enrichi
    filtered_dict = verify_predictions(ref_cds_dict, output_file)           # Vérification des correspondances CDS
    print(f"Nombre de transcripts passés au FASTA = {len(filtered_dict)}")
    analyser_impact(filtered_dict, impact_outfile, transcriptome)     # Analyse de l'impact des mutations
    write_fasta_sequences(filtered_dict, transcriptome, output_fasta, diagnostique_tsv)
    merge_with_uniprot(output_fasta, uniprot_fasta, final_output_fasta)     # Fusion avec UniProt

if __name__ == "__main__":
    main()



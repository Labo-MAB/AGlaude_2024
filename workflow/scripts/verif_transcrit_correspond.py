from Bio import SeqIO

def compare_fasta(file1, file2):
    # Lecture des séquences des deux fichiers FASTA
    seq1 = list(SeqIO.parse(file1, "fasta"))
    seq2 = list(SeqIO.parse(file2, "fasta"))

    # Vérification que les fichiers contiennent le même nombre de séquences
    if len(seq1) != len(seq2):
        print("Les fichiers contiennent un nombre différent de séquences.")
        return
    
    # Comparaison des séquences une à une
    for record1, record2 in zip(seq1, seq2):
        if record1.id != record2.id:
            print(f"Les identifiants sont différents: {record1.id} vs {record2.id}")
            continue

        # Vérification de la longueur des séquences
        if len(record1.seq) != len(record2.seq):
            print(f"Les séquences ont des longueurs différentes pour {record1.id}.")
        
        # Comparaison des séquences lettre par lettre
        differences = []
        for i, (base1, base2) in enumerate(zip(record1.seq, record2.seq)):
            if base1 != base2:
                differences.append((i + 1, base1, base2))  # i + 1 pour position 1-indexée

        # Affichage des résultats
        if differences:
            print(f"Differences in sequence {record1.id}:")
            for pos, base1, base2 in differences:
                print(f"Position {pos}: {base1} -> {base2}")
        else:
            print(f"Aucune différence trouvée dans la séquence {record1.id}.")

# Chemins des fichiers FASTA avec des barres obliques
file1 = "F:/breast_cancer/workflow/data/references/transcriptome.fa"
file2 = "F:/breast_cancer/workflow/results/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/transcripts.fasta"

# Appel de la fonction de comparaison
compare_fasta(file1, file2)

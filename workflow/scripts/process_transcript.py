import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from itertools import chain, combinations, product
import os


def all_subsets(iterable):
    "Retourne toutes les combinaisons non vides des mutations"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))


def parse_vcf(vcf_path):
    mutations = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, pos, _, ref, alts = fields[:5]
            alts = alts.split(",")  # multi-allelique
            for alt in alts:
                mutations.append({
                    "chromosome": chrom,
                    "position": int(pos),
                    "ref": ref,
                    "alt": alt
                })
    return pd.DataFrame(mutations)


def build_transcript_seq_and_map(exons, genome_seq, transcript_id):
    strand = exons.iloc[0]['strand']
    sorted_exons = exons.sort_values(by='start', ascending=(strand == '+'))

    transcript_seq = str(genome_seq.seq)
    genomic_to_transcript = {}
    reconstructed = []
    cumulative = 0

    for _, exon in sorted_exons.iterrows():
        exon_len = exon['end'] - exon['start'] + 1
        if strand == '+':
            for i, pos in enumerate(range(exon['start'], exon['end'] + 1)):
                genomic_to_transcript[pos] = cumulative + i
        else:
            for i, pos in enumerate(range(exon['end'], exon['start'] - 1, -1)):
                genomic_to_transcript[pos] = cumulative + i
        reconstructed.append(transcript_seq[cumulative:cumulative + exon_len])
        cumulative += exon_len

    return "".join(reconstructed), genomic_to_transcript


def apply_mutations(seq, mutations, genomic_to_transcript):
    mutated_seqs = []
    for combo in all_subsets(mutations):
        combo = list(combo)
        valid = True
        for m in combo:
            if m["position"] not in genomic_to_transcript:
                valid = False
                break
            m["transcript_pos"] = genomic_to_transcript[m["position"]]
        if not valid:
            continue

        combo.sort(key=lambda x: x["transcript_pos"])
        new_seq = seq
        offset = 0
        for m in combo:
            pos = m["transcript_pos"] + offset
            ref = m["ref"]
            alt = m["alt"]
            if new_seq[pos:pos + len(ref)] != ref:
                # ref mismatch, skip
                valid = False
                break
            new_seq = new_seq[:pos] + alt + new_seq[pos + len(ref):]
            offset += len(alt) - len(ref)

        if valid:
            signature = "_".join(f"{m['ref']}({m['alt']})" for m in combo)
            mutated_seqs.append((signature, new_seq))
    return mutated_seqs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    genome_seq = next(SeqIO.parse(args.fasta, "fasta"))
    transcript_id = genome_seq.id.split('.')[0]
    exons = pd.read_parquet(args.gtf)
    vcf_df = parse_vcf(args.vcf)

    relevant_mut = vcf_df[vcf_df["chromosome"] == exons.iloc[0]["chromosome"]].copy()
    relevant_mut = relevant_mut[
        (relevant_mut["position"] >= exons["start"].min()) &
        (relevant_mut["position"] <= exons["end"].max())
    ]
    relevant_mut = relevant_mut.to_dict("records")

    transcript_seq, genomic_to_transcript = build_transcript_seq_and_map(exons, genome_seq, transcript_id)
    mutated = apply_mutations(transcript_seq, relevant_mut, genomic_to_transcript)

    records = []
    for signature, new_seq in mutated:
        record_id = f"{transcript_id}_mut_({signature})"
        records.append(SeqRecord(Seq(new_seq), id=record_id, description="Mutated variant"))

    SeqIO.write(records, args.output, "fasta")


if __name__ == "__main__":
    main()

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import argparse


def read_vcf_with_headers(vcf_path):
    """Lit un fichier VCF avec headers et retourne un DataFrame pandas filtré."""
    with open(vcf_path, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]

    # La première ligne non '##' est l'en-tête
    vcf_df = pd.read_csv(
        pd.compat.StringIO("".join(lines)),
        sep="\t"
    )
    vcf_df.rename(columns={vcf_df.columns[0]: "chromosome", vcf_df.columns[1]: "position"}, inplace=True)
    return vcf_df


def extract_transcript_chunks(vcf_path, exons_path, output_dir, n_chunks=32):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Lire le VCF et les exons
    vcf_df = read_vcf_with_headers(vcf_path)
    exons_df = pd.read_parquet(exons_path)

    # Trouver les transcripts impactés par les mutations
    transcript_ids = set()
    for _, row in vcf_df.iterrows():
        chrom = str(row["chromosome"])
        pos = int(row["position"])
        hits = exons_df[
            (exons_df["chromosome"] == chrom) &
            (exons_df["start"] <= pos) &
            (exons_df["end"] >= pos)
        ]
        transcript_ids.update(hits["transcript_id"].unique())

    transcript_ids = sorted(transcript_ids)
    chunks = np.array_split(transcript_ids, n_chunks)

    for i, chunk in enumerate(chunks):
        chunk_file = output_dir / f"chunk_{i:02d}.txt"
        with open(chunk_file, "w") as f:
            for tid in chunk:
                f.write(f"{tid}\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract transcript chunks from VCF and exon data.")
    parser.add_argument("vcf", help="Path to VCF file (with headers)")
    parser.add_argument("exons", help="Path to exons parquet file")
    parser.add_argument("output", help="Directory to write chunk_XX.txt files")
    parser.add_argument("--n_chunks", type=int, default=32, help="Number of chunks to split into")

    args = parser.parse_args()

    extract_transcript_chunks(args.vcf, args.exons, args.output, args.n_chunks)

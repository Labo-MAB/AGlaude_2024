import pandas as pd
import numpy as np
import pyarrow

def extract_exons(gtf_file, output_parquet):
    exons = []
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'exon':
                chromosome = fields[0]
                start, end = int(fields[3]), int(fields[4])
                strand = fields[6]
                attributes = {
                    kv.split(" ")[0]: kv.split(" ")[1].strip('";')
                    for kv in fields[8].split("; ") if kv
                }
                gene_name = attributes.get("gene_name", attributes.get("gene_id", None))
                transcript_id = attributes.get("transcript_id", None)
                gene_biotype = attributes.get("gene_biotype", None)  
                exons.append((chromosome, start, end, strand, gene_name, transcript_id, gene_biotype))
    exons_df = pd.DataFrame(
        exons,
        columns=['chromosome', 'start', 'end', 'strand', 'gene_name', 'transcript_id', 'gene_biotype']
    )
    exons_df[['start', 'end']] = exons_df[['start', 'end']].astype(np.uint32)# Optimisation 
    exons_df['exon_interval'] = pd.IntervalIndex.from_arrays(exons_df['start'], exons_df['end'], closed="both")
    exons_df.to_parquet(output_parquet, index=False)


def main():
    gtf_file = snakemake.input.gtf  
    output_parquet = snakemake.output.parquet
    extract_exons(gtf_file, output_parquet)

if __name__ == "__main__":
    main()

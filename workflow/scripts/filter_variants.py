import sys

def filter_variants(vcf_file, output_file, min_quality=20):
    with open(vcf_file) as vcf, open(output_file, 'w') as out:
        for line in vcf:
            if line.startswith('#'):
                out.write(line)
                continue
            fields = line.strip().split('\t')
            # Quality is typically in the 6th column
            quality = float(fields[5])
            if quality >= min_quality:
                out.write(line)

if __name__ == "__main__":
    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    filter_variants(input_vcf, output_vcf)

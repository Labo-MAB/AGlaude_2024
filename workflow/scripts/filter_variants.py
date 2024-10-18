def filter_variants(vcf_file, output_file, min_quality=20):
    try:
        output_dir = os.path.dirname(output_file)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with open(vcf_file) as vcf, open(output_file, 'w') as out:
            for line in vcf:
                if line.startswith('#'):
                    out.write(line)
                    continue
                fields = line.strip().split('\t')
                try:
                    quality = float(fields[5])
                    if quality >= min_quality:
                        out.write(line)
                except (IndexError, ValueError) as e:
                    print(f"Error processing line: {line.strip()}")
                    print(f"Exception: {e}")
                    raise

    except FileNotFoundError:
        print(f"File not found: {vcf_file}")
        raise
    except IOError as e:
        print(f"I/O error: {e}")
        raise

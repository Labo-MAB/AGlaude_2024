# Exemple de Snakefile pour RNA-Seq

# Réglage des chemins pour les outils et les fichiers
TRIM_GALORE = "trim_galore"
STAR = "STAR"
KALLISTO = "kallisto"
FREEBAYES = "freebayes"
OPENVAR = "OpenVar"

# Réglage des fichiers de référence
REFERENCE_GENOME = "reference/GRCh38.p12.fasta"
TRANSCRIPT_DB = "reference/transcripts.fasta"

# Règle pour le trimming des lectures
rule trim_reads:
    input:
        fq1 = "data/fastq/{sample}_1.fastq.gz",
        fq2 = "data/fastq/{sample}_2.fastq.gz"
    output:
        fq1_trimmed = "data/trimmed/{sample}_1_trimmed.fq.gz",
        fq2_trimmed = "data/trimmed/{sample}_2_trimmed.fq.gz"
    params:
        trim_galore = TRIM_GALORE
    shell:
        "{params.trim_galore} --paired {input.fq1} {input.fq2} --output_dir data/trimmed"

# Règle pour l'alignement avec STAR
rule align_reads:
    input:
        fq1 = "data/trimmed/{sample}_1_trimmed.fq.gz",
        fq2 = "data/trimmed/{sample}_2_trimmed.fq.gz"
    output:
        bam = "data/aligned/{sample}.bam"
    params:
        star = STAR,
        genome = REFERENCE_GENOME
    shell:
        "{params.star} --runThreadN 8 --genomeDir {params.genome} "
        "--readFilesIn {input.fq1} {input.fq2} "
        "--outFileNamePrefix data/aligned/{wildcards.sample} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFilterMismatchNmax 5 --alignSJoverhangMin 10 "
        "--alignMatesGapMax 200000 --alignIntronMax 200000 "
        "--alignSJstitchMismatchNmax '5-1 5 5' "
        "--outSAMprimaryFlag AllBestScore"

# Règle pour la quantification des transcrits avec Kallisto
rule quantify_transcripts:
    input:
        fq1 = "data/trimmed/{sample}_1_trimmed.fq.gz",
        fq2 = "data/trimmed/{sample}_2_trimmed.fq.gz",
        index = "reference/kallisto_index"
    output:
        tpm = "data/quantification/{sample}.tsv"
    params:
        kallisto = KALLISTO
    shell:
        "{params.kallisto} quant --index {input.index} --output-dir data/quantification/{wildcards.sample} "
        "--pairs {input.fq1} {input.fq2}"

# Règle pour l'appel de variants avec FreeBayes
rule call_variants:
    input:
        bam = "data/aligned/{sample}.bam"
    output:
        vcf = "data/variants/{sample}.vcf"
    params:
        freebayes = FREEBAYES
    shell:
        "{params.freebayes} -f {REFERENCE_GENOME} {input.bam} > {output.vcf}"

# Règle pour l'annotation des variants avec OpenVar
rule annotate_variants:
    input:
        vcf = "data/variants/{sample}.vcf",
        transcript_db = TRANSCRIPT_DB
    output:
        annotated_vcf = "data/variants/{sample}_annotated.vcf"
    params:
        openvar = OPENVAR
    shell:
        "{params.openvar} -i {input.vcf} -t {input.transcript_db} -o {output.annotated_vcf}"

# Règle principale
rule all:
    input:
        expand("data/quantification/{sample}.tsv", sample=SAMPLES),
        expand("data/variants/{sample}_annotated.vcf", sample=SAMPLES)

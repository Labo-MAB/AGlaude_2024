# Version 0.0.1 Projet detection variants / nv transcrits du cancer du seins Projet Dre. Brunet
# >>> Ajout  : procedures de manipulation ARNseq data avec les rules
# >>> À ajouter : tout ce qui est path et le bon fonctionnement :)
# Chaque rule est a corriger / rule all a ajouter et rule wildcards 


#merge bam samtool/ module java
###USAGE example(explication/détails des outils): 

#####
#
#
#
#####


configfile: "config.yaml"


import pandas as pd

# Lire les identifiants depuis le fichier texte
with open("SRR_id.txt") as f:
    id_list = f.read().splitlines()

# Lire le fichier TSV pour des comparaisons (exemple d'utilisation, si nécessaire)
df = pd.read_csv('data/comparisons.tsv', sep='\t')

# Définir la règle all pour inclure tous les fichiers de sortie attendus
rule all:
    input:
        expand("data/fastp/{id}_1_trimmed.fastq.gz", id=id_list),
        expand("data/fastp/{id}_2_trimmed.fastq.gz", id=id_list),
        expand("data/aligned/{id}.bam", id=id_list),
        expand("data/quantification/{id}.tsv", id=id_list),
        expand("data/variants/{id}_annotated.vcf", id=id_list)

# Chemins des outils/fichiers
TRIM_GALORE = "trim_galore"
STAR = "STAR"
KALLISTO = "kallisto"
FREEBAYES = "freebayes"
OPENVAR = "OpenVar"

# Fichiers de référence
REFERENCE_GENOME = "reference/GRCh38.p12.fasta"
TRANSCRIPT_DB = "reference/transcripts.fasta"
KALLISTO_INDEX = "reference/kallisto_index"

# Contrôle de qualité FastQC
rule fastqc:
    """Assess the FASTQ quality using FastQC BEFORE TRIMMING"""
    input:
        fq1 = "data/raw_reads/{id}_1.fastq.gz",
        fq2 = "data/raw_reads/{id}_2.fastq.gz"
    output:
        qc_fq1_out = "data/qc/{id}_1_fastqc.html",
        qc_fq2_out = "data/qc/{id}_2_fastqc.html"
    params:
        out_dir = "data/qc"
    log:
        "logs/fastqc_{id}.log"
    threads: 8
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc --outdir {params.out_dir} --format fastq --threads {threads} {input.fq1} {input.fq2} &> {log}"

# Règle pour le trimming
rule trim_reads:
    input:
        fq1 = "data/raw_reads/{id}_1.fastq.gz",
        fq2 = "data/raw_reads/{id}_2.fastq.gz"
    output:
        fq1_trimmed = "data/fastp/{id}_1_trimmed.fastq.gz",
        fq2_trimmed = "data/fastp/{id}_2_trimmed.fastq.gz"
    params:
        trim_galore = TRIM_GALORE
    log:
        "logs/trim_{id}.log"
    shell:
        "{params.trim_galore} --paired {input.fq1} {input.fq2} --output_dir data/fastp --gzip &> {log}"

# Règle pour l'alignement avec STAR
rule align_reads:
    input:
        fq1 = "data/fastp/{id}_1_trimmed.fastq.gz",
        fq2 = "data/fastp/{id}_2_trimmed.fastq.gz"
    output:
        bam = "data/aligned/{id}.bam"
    params:
        star = STAR,
        genome = REFERENCE_GENOME
    log:
        "logs/star_{id}.log"
    threads: 8
    shell:
        "{params.star} --runThreadN {threads} --genomeDir {params.genome} "
        "--readFilesIn {input.fq1} {input.fq2} "
        "--outFileNamePrefix data/aligned/{wildcards.id} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFilterMismatchNmax 5 --alignSJoverhangMin 10 "
        "--alignMatesGapMax 200000 --alignIntronMax 200000 "
        "--alignSJstitchMismatchNmax '5-1 5 5' "
        "--outSAMprimaryFlag AllBestScore &> {log}"

# Règle pour la quantification des transcrits avec Kallisto
rule quantify_transcripts:
    input:
        fq1 = "data/fastp/{id}_1_trimmed.fastq.gz",
        fq2 = "data/fastp/{id}_2_trimmed.fastq.gz",
        index = KALLISTO_INDEX
    output:
        tpm = "data/quantification/{id}.tsv"
    params:
        kallisto = KALLISTO
    log:
        "logs/kallisto_{id}.log"
    shell:
        "{params.kallisto} quant --index {input.index} --output-dir data/quantification/{wildcards.id} "
        "--pairs {input.fq1} {input.fq2} &> {log}"

# Règle pour l'appel de variants avec FreeBayes
rule call_variants:
    input:
        bam = "data/aligned/{id}.bam"
    output:
        vcf = "data/variants/{id}.vcf"
    params:
        freebayes = FREEBAYES
    log:
        "logs/freebayes_{id}.log"
    shell:
        "{params.freebayes} -f {REFERENCE_GENOME} {input.bam} > {output.vcf} 2> {log}"

# Règle pour l'annotation des variants avec OpenVar
rule annotate_variants:
    input:
        vcf = "data/variants/{id}.vcf",
        transcript_db = TRANSCRIPT_DB
    output:
        annotated_vcf = "data/variants/{id}_annotated.vcf"
    params:
        openvar = OPENVAR
    log:
        "logs/openvar_{id}.log"
    shell:
        "{params.openvar} -i {input.vcf} -t {input.transcript_db} -o {output.annotated_vcf} &> {log}"

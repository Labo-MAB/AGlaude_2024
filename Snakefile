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
import os
# Fichiers de référence
REFERENCE_GENOME = "reference/GRCh38.p12.fasta"
TRANSCRIPT_DB = "reference/transcripts.fasta"
KALLISTO_INDEX = "reference/kallisto_index"









# test de paths
fastq_dir = "C:/Users/Anthony/OneDrive - USherbrooke/Documents/1/"
files = [
    "171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L001_1.fastq.gz",
    "171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L001_2.fastq.gz",
    "171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L001_R1_001.220405.A00516.AHVHTNDSX2.fastq.gz"
]

for file in files:
    path = os.path.join(fastq_dir, file)
    print(f"Checking: {path} -> {os.path.exists(path)}")





# Chemins des outils/fichiers
TRIM_GALORE = "trim_galore"
STAR = "STAR"
KALLISTO = "kallisto"
FREEBAYES = "freebayes"
OPENVAR = "OpenVar"

# Définir le répertoire des FASTQ
fastq_dir = "C:/Users/Anthony/OneDrive - USherbrooke/Documents/1/"
output_dir = "data/fastp/"
log_dir = "logs/fastp/"

id_list = [
    "171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L001"
    # Ajoutez d'autres identifiants ici si nécessaire
]

# Assurez-vous que les répertoires de sortie et de logs existent
os.makedirs(output_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

rule all:
    input:
        expand("data/qc/{id}_1_fastqc.html", id=id_list),
        expand("data/qc/{id}_2_fastqc.html", id=id_list),
        #expand("data/fastp/{id}_1_trimmed.fastq.gz", id=id_list),
        #expand("data/fastp/{id}_2_trimmed.fastq.gz", id=id_list),
        expand("data/qc_fastp/{id}_1_fastqc.html", id=id_list),
        expand("data/qc_fastp/{id}_2_fastqc.html", id=id_list),
        #expand("data/aligned/{id}.bam", id=id_list),
        #expand("data/quantification/{id}.tsv", id=id_list),
        #expand("data/variants/{id}_annotated.vcf", id=id_list)

rule fastqc:
    """Assess the FASTQ quality using FastQC BEFORE TRIMMING"""
    input:
        fq1 = fastq_dir + "{id}_R1_001.220405.A00516.AHVHTNDSX2.fastq.gz",
        fq2 = fastq_dir + "{id}_R2_001.220405.A00516.AHVHTNDSX2.fastq.gz"
    output:
        qc_fq1_out = "data/qc/{id}_1_fastqc.html",
        qc_fq2_out = "data/qc/{id}_2_fastqc.html"
    params:
        out_dir = "data/qc"
    log:
        "logs/fastqc_{id}.log"
    threads: 8
    conda:
        "envs/fastqc.yaml"  # a corriger ; conda n'installe pas fastqc
    shell:
        "fastqc --outdir {params.out_dir} --format fastq --threads {threads} {input.fq1} {input.fq2} &> {log}"














# Règle pour le trimming avec fastp
rule trim_reads:
    input:
        fq1 = fastq_dir + "{id}_R1_001.220405.A00516.AHVHTNDSX2.fastq.gz",
        fq2 = fastq_dir + "{id}_R2_001.220405.A00516.AHVHTNDSX2.fastq.gz"
    output:
        fastq1 = "data/fastp/{id}}_R1_001.220405.A00516.AHVHTNDSX2.fastq.gz",
        fastq2 = "data/fastp/{id}_R2_001.220405.A00516.AHVHTNDSX2.fastq.gz",
        unpaired_fastq1 = "data/fastp/{id}_1.unpaired.fastq.gz",
        unpaired_fastq2 = "data/fastp/{id}_2.unpaired.fastq.gz",
        html_report = "data/fastp/{id}_fastp.html",
        json_report = "data/fastp/{id}_fastp.json"

    threads:
        8
    params:
        options = ["--qualified_quality_phred 30", "--length_required 20",
                "--cut_window_size 1", "--cut_mean_quality 30", "--cut_front",
                "--cut_tail"]
    log:
        "logs/fastp/{id}.log"
    conda:
        "../envs/fastp.yml"
    shell:
        "fastp -i {input.fq1} -I {input.fq2} "
        "-o {output.fastq1} -O {output.fastq2} "
        "--unpaired1 {output.unpaired_fastq1} "
        "--unpaired2 {output.unpaired_fastq2} "
        "--thread {threads} "
        "-h {output.html_report} "
        "-j {output.json_report} "
        "{params.options} "
        "&> {log}"


#        
#    params:
#        trim_galore = TRIM_GALORE
#    log:
#        "logs/trim_{id}.log"
#    shell:
#        "{params.trim_galore} --paired {input.fq1} {input.fq2} --output_dir data/fastp --gzip &> {log}"
#
rule qc_fastp:
    """ Assess the FASTQ quality using FastQC AFTER TRIMMING"""
    input:
        trimm_fq1 = rules.trim_reads.output.fastq1,
        trimm_fq2 = rules.trim_reads.output.fastq2,
        trimm_unpaired_fq1 = rules.trim_reads.output.unpaired_fastq1,
        trimm_unpaired_fq2 = rules.trim_reads.output.unpaired_fastq2
    output:
        qc_trimm_fq1_out = "data/qc_fastp/{id}_1_fastqc.html",
        qc_trimm_fq2_out = "data/qc_fastp/{id}_2_fastqc.html"
    params:
        out_dir = "data/qc_fastp"
    log:
        "logs/qc_fastp/{id}.log"
    threads:
        20
    conda:
        "../envs/fastqc.yml"
    shell:
        "fastqc "
        "--outdir {params.out_dir} "
        "--format fastq "
        "--threads {threads} "
        "{input.trimm_fq1} {input.trimm_fq2} "
        "{input.trimm_unpaired_fq1} {input.trimm_unpaired_fq2} "
        "&> {log}" 



# Règle pour l'alignement avec STAR
#rule align_reads:
#    input:
#        fq1 = "data/fastp/{id}_1_trimmed.fastq.gz",
#        fq2 = "data/fastp/{id}_2_trimmed.fastq.gz"
#    output:
#        bam = "data/aligned/{id}.bam"
#    params:
#        star = STAR,
#        genome = REFERENCE_GENOME
#    log:
#        "logs/star_{id}.log"
#    threads: 8
#    shell:
#        "{params.star} --runThreadN {threads} --genomeDir {params.genome} "
#        "--readFilesIn {input.fq1} {input.fq2} "
#        "--outFileNamePrefix data/aligned/{wildcards.id} "
#        "--outSAMtype BAM SortedByCoordinate "
#        "--outFilterMismatchNmax 5 --alignSJoverhangMin 10 "
#        "--alignMatesGapMax 200000 --alignIntronMax 200000 "
#        "--alignSJstitchMismatchNmax '5-1 5 5' "
#        "--outSAMprimaryFlag AllBestScore &> {log}"
#
## Règle pour la quantification des transcrits avec Kallisto
#rule quantify_transcripts:
#    input:
#        fq1 = "data/fastp/{id}_1_trimmed.fastq.gz",
#        fq2 = "data/fastp/{id}_2_trimmed.fastq.gz",
#        index = KALLISTO_INDEX
#    output:
#        tpm = "data/quantification/{id}.tsv"
#    params:
#        kallisto = KALLISTO
#    log:
#        "logs/kallisto_{id}.log"
#    shell:
#        "{params.kallisto} quant --index {input.index} --output-dir data/quantification/{wildcards.id} "
#        "--pairs {input.fq1} {input.fq2} &> {log}"
#
## Règle pour l'appel de variants avec FreeBayes
#rule call_variants:
#    input:
#        bam = "data/aligned/{id}.bam"
#    output:
#        vcf = "data/variants/{id}.vcf"
#    params:
#        freebayes = FREEBAYES
#    log:
#        "logs/freebayes_{id}.log"
#    shell:
#        "{params.freebayes} -f {REFERENCE_GENOME} {input.bam} > {output.vcf} 2> {log}"
#
## Règle pour l'annotation des variants avec OpenVar
#rule annotate_variants:
#    input:
#        vcf = "data/variants/{id}.vcf",
#        transcript_db = TRANSCRIPT_DB
#    output:
#        annotated_vcf = "data/variants/{id}_annotated.vcf"
#    params:
#        openvar = OPENVAR
#    log:
#        "logs/openvar_{id}.log"
#    shell:
#        "{params.openvar} -i {input.vcf} -t {input.transcript_db} -o {output.annotated_vcf} &> {log}"
#
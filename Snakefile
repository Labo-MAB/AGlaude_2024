# Version 0.0.1 Projet detection variants / nv transcrits du cancer du seins Projet Dre. Brunet
# >>> Ajout  : procedures de manipulation ARNseq data avec les rules
# >>> À ajouter : tout ce qui est path et le bon fonctionnement :)
# Chaque rule est a corriger / rule all a ajouter et rule wildcards 


#merge bam samtool/ module java
###USAGE example:

#



###


#configfile 


#  des chemins des outils/fichiers
TRIM_GALORE = "trim_galore"
STAR = "STAR"
KALLISTO = "kallisto"
FREEBAYES = "freebayes"
OPENVAR = "OpenVar"

# fichiers de référence 
REFERENCE_GENOME =             #"reference/GRCh38.p12.fasta"     ## l'article customDB
TRANSCRIPT_DB =                #"reference/transcripts.fasta"

# controle de qualité FastQC

rule fastqc:
    """ Assess the FASTQ quality using FastQC BEFORE TRIMMING"""
    input:
        fq1 = "/{id}_1.fastq.gz",
        fq2 = "/{id}_2.fastq.gz"
    output:
        qc_fq1_out = "/qc/{id}_1_fastqc.html",
        qc_fq2_out = "/qc/{id}_2_fastqc.html"
    params:
        out_dir = "/qc"
    log:
        
    threads:
        
    conda:
        
    shell:
        "fastqc "
        "--outdir {params.out_dir} "
        "--format fastq "
        "--threads {threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}" 





# À ajouter Controle de qualité   thread 32?



rule all:
    input:
        expand(),
        expand()


# Règle pour le trimming   trim_galore   # Vérifier l'outil
rule trim_reads:
    input:
        fq1 = "/{id}_1.fastq.gz",
        fq2 = "/{id}_2.fastq.gz"
    output:
        fq1_trimmed = "/{id}_1_trimmed.fq.gz",
        fq2_trimmed = "/{id}_2_trimmed.fq.gz"
    params:
        trim_galore = TRIM_GALORE
    shell:
        "{params.trim_galore} --paired {input.fq1} {input.fq2} --output_dir data/trimmed"



# À ajouter controle de qualité (servira de vérification)
 
# Règle pour l'alignement  STAR     # Vérifier l'outil
rule align_reads:
    input:
        fq1 = "/{id}_1_trimmed.fq.gz",
        fq2 = "/{id}_2_trimmed.fq.gz"
    output:
        bam = "data/aligned/{id}.bam"
    params:
        star = STAR,
        genome = REFERENCE_GENOME
    shell:
        "{params.star} --runThreadN 8 --genomeDir {params.genome} "
        "--readFilesIn {input.fq1} {input.fq2} "
        "--outFileNamePrefix data/aligned/{wildcards.id} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFilterMismatchNmax 5 --alignSJoverhangMin 10 "
        "--alignMatesGapMax 200000 --alignIntronMax 200000 "
        "--alignSJstitchMismatchNmax '5-1 5 5' "
        "--outSAMprimaryFlag AllBestScore"


#--runMode genomeGenerate
#--genomeDir /path/to/genomeDir
#--genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ...
#--sjdbGTFfile /path/to/annotations.gtf
#--sjdbOverhang ReadLength-1

# Règle pour la quantification des transcrits avec Kallisto    # Vérifier l'outil
rule quantify_transcripts:
    input:
        fq1 = "data/trimmed/{id}_1_trimmed.fq.gz",
        fq2 = "data/trimmed/{id}_2_trimmed.fq.gz",
        index = "reference/kallisto_index"
    output:
        tpm = "data/quantification/{id}.tsv"
    params:
        kallisto = KALLISTO
    shell:
        "{params.kallisto} quant --index {input.index} --output-dir data/quantification/{wildcards.id} "
        "--pairs {input.fq1} {input.fq2}"

# Règle pour l'appel de variants avec FreeBayes     # Vérifier l'outil
rule call_variants:
    input:
        bam = "data/aligned/{id}.bam"
    output:
        vcf = "data/variants/{id}.vcf"
    params:
        freebayes = FREEBAYES
    shell:
        "{params.freebayes} -f {REFERENCE_GENOME} {input.bam} > {output.vcf}"

# Règle pour l'annotation des variants avec OpenVar    # Vérifier l'outil
rule annotate_variants:
    input:
        vcf = "data/variants/{id}.vcf",
        transcript_db = TRANSCRIPT_DB
    output:
        annotated_vcf = "data/variants/{id}_annotated.vcf"
    params:
        openvar = OPENVAR
    shell:
        "{params.openvar} -i {input.vcf} -t {input.transcript_db} -o {output.annotated_vcf}"


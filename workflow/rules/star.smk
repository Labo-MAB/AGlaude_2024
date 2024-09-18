rule star_index:
    """ Generates the genome index for STAR """
    input:
        fasta = rules.download_human_genome.output.genome,
        gtf = config['download']['human_gtf']
    output:
        chrNameLength = config['path']['chrNameLength']
    params:
        dir = config['path']['star_index']
    log:
        "logs/STAR/index.log"
    conda:
        "../envs/star.yml"
    threads:
        32
    shell:
        "mkdir -p {params.dir} && "
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99"
        "&> {log}"


rule star_alignReads:
    """ Generates a bam file using STAR """
    input:
        idx = rules.star_index.output,
        fq1 = rules.trim_reads.output.gal_trim1,
        fq2 = rules.trim_reads.output.gal_trim2
    output:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam",
        bam_logs = "results/STAR/{id}/Log.final.out"
    params:
        index = config['path']['star_index'],
        output_dir = "results/STAR/{id}/"
    log:
        "logs/STAR/{id}.log"
    threads:
        32
    conda:
        "../envs/star.yml"
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.index} "
        "--readFilesIn {input.fq1} {input.fq2}  "
        "--runThreadN {threads} "
        "--readFilesCommand zcat "
        "--outReadsUnmapped Fastx "
        "--outFilterType BySJout "
        "--outStd Log "
        "--outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.output_dir} "
        "--outFilterScoreMinOverLread 0.3 "
        "--outFilterMatchNminOverLread 0.3 "
        "--outFilterMultimapNmax 100 "
        "--winAnchorMultimapNmax 100 "
        "--limitBAMsortRAM 600000000000"
        "--alignEndsProtrude 5 ConcordantPair "
        "&> {log}"


############# Star se fera dans une config a part 
# Règle pour l'alignement avec STAR
#rule align_reads:
#    input:
#        fq1 = "data/trim_galore/{id}_1_trimmed.fastq.gz",
#        fq2 = "data/trim_galore/{id}_2_trimmed.fastq.gz"
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

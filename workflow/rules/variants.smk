rule index_genome:
    input:
        genome=rules.download_human_genome.output.genome
    output:
        fai="data/references/genome_fa/homo_sapiens_genome.fa.fai"
    conda:
        "../envs/freebayes.yml"
    shell:
        "samtools faidx {input.genome}"

rule create_targets:
    input:
        fai=rules.index_genome.output.fai
    output:
        targets="results/variants/{id}/targets.txt"
    conda:
        "../envs/freebayes.yml"
    shell:
        "scripts/fasta_generate_regions.py {input.fai} 1000000 > {output.targets}"

rule index_bam:
    input:
        bam=rules.star_alignreads.output.bam
    output:
        bai="results/STAR/{id}/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/freebayes.yml"
    shell:
        "samtools index {input.bam}"

rule call_variants:
    input:
        bam=rules.star_alignreads.output.bam,
        bai=rules.index_bam.output.bai,
        genome=rules.download_human_genome.output.genome,
        fai=rules.index_genome.output.fai,
        targets=rules.create_targets.output.targets
    output:
        vcf="results/variants/{id}/variants.vcf"
    params:
        min_alternate_count=5,
        min_coverage=10
    conda:
        "../envs/freebayes.yml"
    log:
        "logs/freebayes_{id}/freebayes.log"
    threads: 32  
    shell:
        """
        freebayes-parallel {input.targets} {threads} -f {input.genome} {input.bam} > {output.vcf} 2> {log}
        """

rule filter_variants:
    input:
        vcf=rules.call_variants.output.vcf
    output:
        vcf_filtered="results/variants/{id}/20QC_variant.vcf" 
    conda:
        "../envs/bcftools.yml"  
    log:
        "logs/freebayes_{id}/filter_variants.log"
    shell:
        "bcftools filter -s LowQual -e 'QUAL<20' {input.vcf} -o {output.vcf_filtered} > {log} 2>&1"

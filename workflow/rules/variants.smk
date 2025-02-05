rule call_variants:  # Freebayes pour avoir la liste_mutant.vcf des dans l'ARNseq
    input:
        bam = rules.star_alignreads.output.bam,
        genome = rules.download_human_genome.output.genome,
        gtf = rules.download_human_gtf.output.gtf
    output:
        vcf = "results/variants/{id}/variants.vcf",
        annotated_vcf = "results/variants/{id}/variants_annotated.vcf"
    params:
        out_dir = "results/variants",
        min_alternate_count = 5,
        min_coverage = 10
    conda:
        "../envs/freebayes.yml"
    log:
        "logs/freebayes_{id}/freebayes.log"
    threads: 8  
    shell: 
        """
        set -e  # Arrête l'exécution en cas d'erreur
        mkdir -p results/variants && \
        freebayes -f {input.genome} \
            --min-alternate-count {params.min_alternate_count} \
            --min-coverage {params.min_coverage} \
            {input.bam} \
            > {output.vcf} \
            2>> {log} && \
        snpeff -hgvsTrId -geneId hg38 {output.vcf} -o vcf > {output.annotated_vcf} 2>> {log}
        """

rule filter_variants: # QC des mutants avec trashold de 20 
    input:
        vcf = rules.call_variants.output.annotated_vcf  
    output:
        vcf_filtered = "results/variants/{id}/20QC_variant.vcf" 
    conda:
        "../envs/bcftools.yml"  
    log:
        "logs/freebayes_{id}/filter_variants.log"
    shell:
        "bcftools filter -s LowQual -e 'QUAL<20' {input.vcf} -o {output.vcf_filtered} > {log} 2>&1"


rule apply_variants:
    input:
        fasta = rules.download_human_genome.output.genome,  
        vcf = rules.filter_variants.output.vcf_filtered,  
        gtf = rules.download_human_gtf.output.gtf, 
        transcriptome = rules.build_filtered_transcriptome.output.transcriptome_final_custom
    output:
        fasta = "results/final/{id}/transcrits_variants.fa"  
    conda:
        "../envs/python.yml"  
    log:
        "logs/freebayes_{id}/apply_variants.log" 
    script:
        "../scripts/add_variants.py"  

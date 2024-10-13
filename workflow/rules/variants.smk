rule call_variants:
    input:
        bam = rules.star_alignreads.output.bam,
        genome = rules.download_human_genome.output.genome
    output:
        vcf = "results/variants/{id}.vcf"
    params:
        out_dir = "results/variants",
        min_alternate_count = 0.5,  
        min_coverage = 10
    conda:
        #"/home/anthony/miniconda3/envs/freebayes"
        "../envs/freebaye.yml"
    log:
        "logs/freebayes_{id}.log"
    shell: 
        """
        mkdir -p {params.out_dir} &&
        freebayes -f {input.genome} --min-alternate-fraction {params.min_alternate_count} 
        --min-coverage {params.min_coverage} {input.bam} > {output.vcf} 2> {log}
        """


rule filter_variants:
    input:
        vcf = rules.call_variants.output.vcf
    output:
        vcf_filtered = "data/variants/{id}_filtred.vcf"
    conda:
        "../envs/python.yml"  # Environnement contenant Python
    log:
        "logs/filter_variants_{id}.log"
    script:
        "../scripts/filter_variants.py"  
# Règle pour l'annotation des variants avec OpenVar
#rule annotate_variants:
#    input:
#        vcf = rules.call_variants.output.vcf
#    output:
#        annotated_vcf = "data/variants/{id}_annotated.vcf"
#    conda:
#        "../envs/openvar.yml"  # Assurez-vous que ce chemin est correct
#    log:
#        "logs/openvar_{id}.log"
#    shell:
#        """
#        openvar -i {input.vcf} -o {output.annotated_vcf} &> {log}
#        """

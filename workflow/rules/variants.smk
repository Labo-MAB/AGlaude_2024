rule call_variants: # Freebayes pour avoir la liste_mutant.vcf des dans l'ARNseq
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
        "logs/freebayes_{id}.log"
    threads: 8  
    shell: 
        """
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
        "logs/filter_variants_{id}.log"
    shell:
        #"scripts/filter_variants.py" 
        "bcftools filter -s LowQual -e 'QUAL<20' {input.vcf} -o {output.vcf_filtered} > {log} 2>&1"



# Next, il faut que un fichier transcriptome non mutant avec tout les gènes ayant un mutant 
# input -> transcriptome_custom + filter_variants avec grep 

rule filter_transcriptome_by_mutations:
    input:
        transcriptome = rules.build_filtered_transcriptome.output.transcriptome_final_custom,
        vcf = rules.filter_variants.output.vcf_filtered 
    output:
        filtered_transcriptome = "results/{id}/filtered_transcriptome_with_mutations.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/filter_mutated_genes_and_filter_transcriptome.py"  # Script Python pour extraire les gènes mutés et filtrer directement


#rule apply_variants: 
#    input:
#        #fasta = rules.build_transcriptome.output.transcriptome,
#        fasta = rules.download_human_genome.output.genome,  # Remplacez par le chemin réel
#        vcf = rules.filter_variants.output.vcf_filtered,  # Fichier VCF filtré
#        gtf = rules.download_human_gtf.output.gtf  # Fichier GTF
#    output:
#        fasta = "results/variants/{id}/transcrits_variants.fa"  # Fichier de sortie pour les transcrits avec variantes
#    conda:
#        "../envs/python.yml"  # Environnement Conda si nécessaire
#    log:
#        "logs/apply_variants_{id}.log"
#    script:
#        "../scripts/add_variants.py"  # Chemin vers votre script
#
#
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
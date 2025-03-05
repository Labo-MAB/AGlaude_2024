rule call_variants:
    input:
        bam=rules.star_alignreads.output.bam,
        genome=rules.download_human_genome.output.genome,
        gtf=rules.download_human_gtf.output.gtf
    output:
        vcf="results/variants/{id}/variants.vcf",
    params:
        min_alternate_count=5,
        min_coverage=10,
    conda:
        "../envs/freebayes.yml"
    log:
        "logs/freebayes_{id}/freebayes.log"
    threads: 8  
    shell: 
        "freebayes -f {input.genome} {input.bam} > {output.vcf} 2> {log}"


        #java -jar {input.snpeff} -hgvsTrId -geneId hg38 {output.vcf} -o vcf > {output.annotated_vcf} 2>> {log}  
        #"""
        #set -e  # Arrête l'exécution en cas d'erreur
        ## Appel de variants avec FreeBayes
        #freebayes -f {input.genome} \
#
        #    {input.bam} \
        #    > {output.vcf} 2>> {log}
        #"""
rule VCF_annotation:
    input:
        vcf=rules.call_variants.output.vcf
    output:
        annotated_vcf="results/variants/{id}/variants_annotated.vcf"
    params:
        SNPEFF_JAR="$EBROOTSNPEFF/snpEff.jar",
        SNPEFF_DB="data/references/data/hg38/snpEffectPredictor.bin", #  a corriger
        CONFIG="data/references/snpEff.config" # a corriger
    #conda:
    #    "../envs/freebayes.yml"
    log:
        "logs/freebayes_{id}/annotation.log"
    threads: 8  
    shell: 
        """
        module load snpeff
        module load java
        java -jar {params.SNPEFF_JAR} -c {params.CONFIG} -hgvsTrId -geneId hg38 -noDownload {input.vcf} > {output.annotated_vcf} 2> {log}
        """

        #snpeff -hgvsTrId -geneId hg38 {input.vcf} > {output.annotated_vcf} 2>> {log}


rule filter_variants: # QC des mutants avec trashold de 20 
    input:
        vcf = rules.VCF_annotation.output.annotated_vcf  
    output:
        vcf_filtered = "results/variants/{id}/20QC_variant.vcf" 
    conda:
        "../envs/bcftools.yml"  
    log:
        "logs/freebayes_{id}/filter_variants.log"
    shell:
        "bcftools filter -s LowQual -e 'QUAL<20' {input.vcf} -o {output.vcf_filtered} > {log} 2>&1"

rule build_exon_dataframe:
    """Extracts exons from the GTF file and saves them in a Parquet file."""
    input:
        gtf = rules.download_human_gtf.output.gtf  
    output:
        parquet = "results/exon_data.parquet"  
    conda:
        "../envs/python.yml"  
    script:
        "../scripts/build_exon_dataframe.py"

rule apply_variants:
    input:
        exon_parquet = rules.build_exon_dataframe.output.parquet,
        #fasta = rules.download_human_genome.output.genome,  
        vcf = rules.filter_variants.output.vcf_filtered,  
        transcripts_fasta = rules.build_filtered_transcriptome.output.transcriptome_final_custom
    output:
        mutated_output_fasta = "results/final/{id}/mutated_transcripts.fa",
        combined_output_fasta = "results/final/{id}/combined_transcripts.fa"
    conda:
        "../envs/python.yml"  
    script:
        "../scripts/add_variants.py" 

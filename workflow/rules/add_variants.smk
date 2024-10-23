#rule vep_annotation:
#    input:
#        vcf= rules.call_variants.output.vcf
#    output:
#        vcf="results/variants/{id}/annotated_variants.vcf"
#    conda:
#        "../envs/vep.yml"
#    shell:
#        """
#        vep --input_file {input.vcf} --output_file {output.vcf} --format vcf --vcf --symbol --canonical --hgvs --protein --ccds --uniprot --biotype --nearest symbol --distance 5000 --species homo_sapiens 2> error.log
#        """
rule bgzip_vcf:
    input:
        vcf= rules.filter_variants.output.vcf_filtered,
    output:
        vcf_gz="results/variants/{id}/variant_filtered.vcf.gz",
        vcf_gz_tbi="results/variants/{id}/variant_filtered.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/bgzip_vcf_{id}.log"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        """
rule bcftools_consensus:
    input:
        vcf = rules.bgzip_vcf.output.vcf_gz,
        ref = rules.download_human_genome.output.genome
    output:
        fasta = "results/variants/{id}/mutated_sequence.fasta"
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/bcftools_{id}.log"
    shell:
        """
        bcftools consensus -f {input.ref} {input.vcf} > {output.fasta} 2> {log}
        """


rule merge_txvariant:
    input:
        vcf = rules.call_variants.output.vcf,
        abundance_tsv = "results/dge/kallisto/{id}/abundance.tsv"
    output:
        fasta = "results/merged/{id}/{id}.fasta"
    conda:
        "../envs/add_variants.yml"  
    log:
        "logs/add_variants_{id}.log"
    script:
        "../scripts/add_variants.py"

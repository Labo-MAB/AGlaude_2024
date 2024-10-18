rule vep_annotation:
    input:
        vcf= rules.call_variants.output.vcf
    output:
        vcf="results/variants/{id}/annotated_variants.vcf"
    conda:
        "../envs/vep.yml"
    shell:
        """
        vep --input_file {input.vcf} --output_file {output.vcf} --format vcf --vcf --symbol --canonical --hgvs --protein --ccds --uniprot --biotype --nearest symbol --distance 5000 --species homo_sapiens
        """

rule bcftools_consensus:
    input:
        vcf= rules.vep_annotation.output.vcf,
        ref= rules.download_human_genome.output.genome
    output:
        fasta="results/{id}/mutated_sequence.fasta"
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools consensus -f {input.ref} {input.vcf} > {output.fasta}
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

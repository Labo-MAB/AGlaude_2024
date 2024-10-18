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

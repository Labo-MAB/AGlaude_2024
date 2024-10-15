rule merge_txvariant: 
    input:
        vcf = rules.call_variants.output.vcf
        abundance_tsv = rules.kalisto_quant.output
    output:
        fasta = "results/merged/{id}/{id}_filtred.fasta"
    conda:
        "../envs/python.yml"  # Environnement contenant Python
    log:
        "logs/filter_variants_{id}.log"
    script:
        "../scripts/filter_variants.py"  
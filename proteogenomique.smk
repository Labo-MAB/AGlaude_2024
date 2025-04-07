rule geno2proteo:
    input:
        mutated_transcriptome=rules.filter_transcripts.output.filtered_fasta
    output:
        proteotome_personalised="result/final/{id}/proteotome_personalised.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/bob.py"

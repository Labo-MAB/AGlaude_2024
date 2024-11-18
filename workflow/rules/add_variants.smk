



#rule bgzip_gtf:
#    input:
#        gtf = rules.download_human_gtf.output.gtf
#    output:
#        gtf_gz="data/references/homo_sapiens_{id}.gtf.gz",
#        gtf_idx="data/references/homo_sapiens_{id}.gtf.gz.tbi"  
#    conda:
#        "../envs/bcftools.yml"
#    log:
#        "logs/bgzip_gtf_{id}.log"
#    shell:
#        """
#        bgzip -c {input.gtf} > {output.gtf_gz}
#        tabix -p gff {output.gtf_gz} 
#        """
#
#
#rule vep_annotation:
#    input:
#        vcf = rules.filter_variants.output.vcf_filtered,
#        fasta = rules.download_human_genome.output.genome,
#        gtf_gz = rules.bgzip_gtf.output.gtf_gz,  # Le GTF compressé
#        gtf_idx = rules.bgzip_gtf.output.gtf_idx   # L'index
#    output:
#        vcf="results/variants/{id}/annotat_variants.vcf"
#    conda:
#        "../envs/vep.yml"
#    log:
#        "logs/vep_{id}.log"
#    shell:
#        """
#        vep --input_file {input.vcf} --output_file {output.vcf} --format vcf --vcf \
#        --symbol --canonical --hgvs --protein --ccds --uniprot --biotype \
#        --nearest symbol --distance 5000 --species homo_sapiens \
#        --fasta {input.fasta} --gtf_gz {input.gtf_gz} 2> error.log 
#        """
#
#
#
#rule apply_mutations_to_transcripts:
#    input:
#        vcf=rules.vep_annotation.output.vcf,  # fichier VCF annoté (non compressé)
#        transcripts_fasta=rules.build_transcriptome.output  # transcrits extraits
#    output:
#        mutated_transcripts_fasta="results/variants/{id}/test_transcripts.fa"  # transcrits avec mutations appliquées
#    conda:
#        "../envs/bcftools.yml"  
#    log:
#        "logs/apply_mutations_to_transcripts_{id}.log"
#    shell:
#        """
#        bcftools consensus -f {input.transcripts_fasta} -o {output.mutated_transcripts_fasta} {input.vcf} 2> {log}
#        """
#
#
#
#
#rule bgzip_vcf:
#    input:
#        vcf= rules.filter_variants.output.vcf_filtered,
#    output:
#        vcf_gz="results/variants/{id}/variant_filtered.vcf.gz",
#        vcf_gz_tbi="results/variants/{id}/variant_filtered.vcf.gz.tbi"
#    conda:
#        "../envs/bcftools.yml"
#    log:
#        "logs/bgzip_vcf_{id}.log"
#    shell:
#        """
#        bgzip -c {input.vcf} > {output.vcf_gz}
#        tabix -p vcf {output.vcf_gz}
#        """
##rule bcftools_mutant_genom:
##    input:
##        vcf = rules.bgzip_vcf.output.vcf_gz,
##        ref = rules.download_human_genome.output.genome
#    output:
#        fasta = "results/variants/{id}/mutated_sequence.fasta"
#    conda:
#        "../envs/bcftools.yml"
#    log:
#        "logs/bcftools_{id}.log"
#    shell:
#        """
#        bcftools consensus -f {input.ref} {input.vcf} > {output.fasta} 2> {log}
#        """
#
#
#rule merge_txvariant:
#    input:
#        vcf = rules.call_variants.output.vcf,
#        abundance_tsv = "results/dge/kallisto/{id}/abundance.tsv",
#        ref_fasta = rules.bcftools_mutant_genom.output.fasta,  
#        gtf =  rules.download_human_gtf.output.gtf
#    output:
#        fasta = "results/variants/{id}/transcripts.fasta"  
#    conda:
#        "../envs/gffread.yml"  
#    log:
#        "logs/merge_txvariant_{id}.log"
#    shell:
#        """
#        gffread {input.gtf} -g {input.ref_fasta} -w {output.fasta}
#        """
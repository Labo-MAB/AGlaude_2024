CHR_7_TO_15 = [str(c) for c in range(7, 16)]
CHR_16_TO_22_XYMT = [str(c) for c in range(16, 23)] + ["X", "Y", "MT"]

rule build_exon_dataframe:
    """Extracts exons from the GTF file and saves them in a Parquet file."""
    input:
        gtf=rules.download_human_gtf.output.gtf  
    output:
        parquet="data/references/exon_data.parquet"  
    conda:
        "../envs/python.yml"  
    script:
        "../scripts/build_exon_dataframe.py"

rule split_vcf_1_2:
    input:
        vcf=rules.filter_variants.output.vcf_filtered
    output:
        vcf="results/{id}/split_vcf/1_2.vcf"
    params:
        chromosomes=["1", "2"]
    conda:
        "../envs/bcftools.yml"  
    script:
        "../scripts/split_vcf_by_group.py"

rule split_vcf_3_4_5_6:
    input:
        vcf=rules.filter_variants.output.vcf_filtered
    output:
        vcf="results/{id}/split_vcf/3_4_5_6.vcf"
    params:
        chromosomes=["3", "4", "5", "6"]
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/split_vcf_by_group.py"

rule split_vcf_7_to_15:
    input:
        vcf=rules.filter_variants.output.vcf_filtered
    output:
        vcf="results/{id}/split_vcf/7_to_15.vcf"
    params:
        chromosomes=CHR_7_TO_15
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/split_vcf_by_group.py"

rule split_vcf_rest:
    input:
        vcf=rules.filter_variants.output.vcf_filtered
    output:
        vcf="results/{id}/split_vcf/rest.vcf"
    params:
        chromosomes=CHR_16_TO_22_XYMT,
        include_rest=True
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/split_vcf_by_group.py"

rule apply_variants12:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_1_2.output.vcf,
        transcripts_fasta=rules.build_filtered_transcriptome.output.transcriptome_final_custom
    output:
        mutated_output_fasta="results/final/{id}/12_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule apply_variants3456:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_3_4_5_6.output.vcf,
        transcripts_fasta=rules.build_filtered_transcriptome.output.transcriptome_final_custom
    output:
        mutated_output_fasta="results/final/{id}/3456_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule apply_variants7to15:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_7_to_15.output.vcf,
        transcripts_fasta=rules.build_filtered_transcriptome.output.transcriptome_final_custom
    output:
        mutated_output_fasta="results/final/{id}/7to15_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule apply_variants_rest:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_rest.output.vcf,
        transcripts_fasta=rules.build_filtered_transcriptome.output.transcriptome_final_custom
    output:
        mutated_output_fasta="results/final/{id}/rest_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule fusion_transcripts:
    input:
        mutated_output_fasta1=rules.apply_variants12.output.mutated_output_fasta,
        mutated_output_fasta2=rules.apply_variants3456.output.mutated_output_fasta,
        mutated_output_fasta3=rules.apply_variants7to15.output.mutated_output_fasta,
        mutated_output_fasta4=rules.apply_variants_rest.output.mutated_output_fasta,
        transcripts_fasta=rules.build_filtered_transcriptome.output.transcriptome_final_custom
    output:
        combined_transcripts="results/final/{id}/combined_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/combine_transcripts.py"

rule filter_transcripts:
    "Fichier de sortie pour la portion proteomique"
    input:
        combined_fasta=rules.fusion_transcripts.output.combined_transcripts
    output:
        filtered_fasta="results/final/{id}/filtered_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/filter_transcripts.py"

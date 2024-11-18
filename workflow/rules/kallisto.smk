# Règle pour la quantification des transcrits avec Kallisto  **** A corriger, les paths && les dependances entre chaque rules)

rule build_transcriptome:  # Sert de transcriptome de reference avec tout les id transcrits
    input:
        genome = rules.download_human_genome.output.genome,
        gtf = rules.download_human_gtf.output.gtf
    output:
        transcriptome = config["path"]["transcriptome"]
    conda:
        "../envs/gffread.yml"
    message:
        "Build a reference transcriptome using gffread."
    log:
        "logs/build_transcriptome/build_transcriptome.log"
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output}"  


rule kallisto_index: # L'index qui sert a quantifier
    input:
        rules.build_transcriptome.output
    output:
        "data/references/kallisto.idx"
    conda:
        "../envs/kallisto.yml"
    log:
        "logs/kallisto/index.log"
    message:
        "Builds an index from the FASTA file."
    shell:
        "kallisto index "
        "--index={output} "
        "{input} "
        "&> {log}"



rule kallisto_quant: # Donne labandonce des transcrits chez l'individu
    input:
        idx = rules.kallisto_index.output,
        fq1 = rules.trim_reads.output.gal_trim1,
        fq2 = rules.trim_reads.output.gal_trim2
    output: 
        abundance_tsv = "results/dge/kallisto/{id}/abundance.tsv"
    params:
        bootstrap = "50",
        outdir = "results/dge/kallisto/{id}/"
    threads:
        1
    conda:
        "../envs/kallisto.yml"
    log:
        "logs/kallisto/{id}.log"
    message:
        "Perform pseudoalignment and quantify transcript abundance for {wildcards.id}."
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}"


rule filter_abundance: # ne contient que les transcrits présents dans l'ARNseq (donc retire les 0) Ajouter un trashold pour le nombre de lecture pour etre accepter
    input:
        abundance = rules.kallisto_quant.output.abundance_tsv
    output:
        filtered_abundance = "results/dge/kallisto/{id}/abundance_filtered.tsv",  # Contient seulement les gènes quantifiers 
        transcript_ids = "results/dge/kallisto/{id}/transcript_ids.txt" # contient seulement le ID des gènes quantifiers
    log:
        "logs/kallisto/filter_{id}.log"
    shell:
        "awk '$4 > 0' {input.abundance} > {output.filtered_abundance} && "
        "cut -f1 {output.filtered_abundance} > {output.transcript_ids}"


rule build_filtered_transcriptome:
    input:
        gtf = rules.download_human_gtf.output.gtf,
        genome = rules.download_human_genome.output.genome,
        transcript_ids = rules.filter_abundance.output.transcript_ids
    output:
        filtered_gtf = "results/{id}/filtered_transcripts.gtf",
        transcriptome_final_custom = "results/{id}/transcriptome_final_custom.fa"
    conda:
        "../envs/gffread.yml"
    message:
        "Build a filtered reference transcriptome using gffread."
    log:
        "logs/build_transcriptome/build_filtered_transcriptome_{id}.log"
    shell:
        """
        # Filtrer le GTF pour ne garder que les transcrits désirés
        grep -F -f {input.transcript_ids} {input.gtf} > {output.filtered_gtf}
        
        # Modifier le GTF pour remplacer 'geneID' par 'gene_name'
        awk '{{ if ($0 ~ /geneID/) {{ gsub("geneID", "gene_name"); }} print $0; }}' {output.filtered_gtf} > {output.filtered_gtf}_modified.gtf
        
        # Générer le transcriptome à partir du GTF filtré et modifié
        gffread {output.filtered_gtf}_modified.gtf -g {input.genome} -w {output.transcriptome_final_custom}
        """


# À Utiliser la 1st fois, ensuite le fichier est sauvegardé
rule download_human_gtf:
    """ Download gtf of human genome from Ensembl """
    output:
        gtf = 'data/references/gtf/homo_sapiens.gtf'
    params:
        link = config['download']['human_gtf']
    shell:
        #"mkdir -p data/references/gtf &&"
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"



rule download_human_gff3:
    """ Download gff3 of human genome from Ensembl """
    output:
        gff3 = 'data/references/gff3/homo_sapiens.gff3'
    params:
        link = config['download']['human_gff3']
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gff3}"


rule download_human_genome:
    """Download the reference genome (fasta file) of human
       from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/homo_sapiens_genome.fa'
    params:
        link = config['download']['human_genome_fa']
    shell:
        "mkdir -p data/references/genome_fa && " 
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

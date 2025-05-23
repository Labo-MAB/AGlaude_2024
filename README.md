- Éxecution du pipeline Snakemake
Le pipeline doit être lancé à partir du répertoire workflow dans l'environnement conda où vous avez installé snakemake.

Dans les récentes versions de Snakemake (à partir de 7), vous pouvez vous créer 2 types de profil pour lancer le pipeline:

profile_local: Si les nœuds du cluster n'ont pas accès à Internet, exécutez en premier lieu les tâches nécessitant Internet localement (i.e all downloads) depuis le répertoire du workflow avec :
# Lancer la section download_genome du pipeline localement
snakemake download_genome --profile ../profile_local/
# Lancer un workflow localement
snakemake --profile ../profile_local/ 

# Lancer un workflow sur un cluster (ou section all du Snakefile)
snakemake --profile ../profile_slurm/

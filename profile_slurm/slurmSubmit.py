#!/usr/bin/env python3

"""
Adapted from https://bitbucket.org/snakemake/snakemake/issues/28/clustering-jobs-with-snakemake

Launch with :
snakemake -j 99 --use-conda --cluster-config cluster.json --immediate-submit --notemp --cluster 'python3 slurmSubmit.py {dependencies}'

"""
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

cmdline = "sbatch "
for param, val in job_properties['cluster'].items():
    cmdline += "--{param} {val} ".format(param=param, val=val)

# Set up dependencies
dependencies = set(sys.argv[1:-1])
if dependencies:
    cmdline += " --dependency=afterok:{} ".format(":".join(dependencies))

# Adding the actual job
cmdline += jobscript

# remove the leading and trailing white space for the submitted jobid
cmdline += r" | awk '{print substr($NF, 0, length($NF))}'"

sys.stdout.write(cmdline)

os.system(cmdline)

#!/usr/bin/env python3
#
#import os
#import sys
#import subprocess
#from snakemake.utils import read_job_properties
#
## Lire le script de job généré par Snakemake
#jobscript = sys.argv[-1]
#job_properties = read_job_properties(jobscript)
#
## Construction de la commande sbatch
#cmdline = ["sbatch"]
#
## Ajout des options de ressources si définies dans le Snakefile
#cluster_params = job_properties.get("cluster", {})
#for param, val in cluster_params.items():
#    cmdline.append(f"--{param}={val}")
#
## Gestion des dépendances entre les jobs
#dependencies = sys.argv[1:-1]
#if dependencies:
#    dependency_str = ":".join(dependencies)
#    cmdline.append(f"--dependency=afterok:{dependency_str}")
#
## Ajout du script de job à la fin de la commande
#cmdline.append(jobscript)
#
## Soumission du job avec subprocess pour capturer proprement l'ID du job SLURM
#try:
#    result = subprocess.run(cmdline, stdout=subprocess.PIPE, check=True, text=True)
#    # Extraire l'ID du job SLURM (normalement le dernier mot de la sortie de sbatch)
#    job_id = result.stdout.strip().split()[-1]
#    print(job_id)
#except subprocess.CalledProcessError as e:
#    print(f"Erreur lors de la soumission du job : {e}", file=sys.stderr)
#    sys.exit(1)
# am pm
#
# #!/usr/bin/env python3

# """
# Adapted from https://bitbucket.org/snakemake/snakemake/issues/28/clustering-jobs-with-snakemake

# Launch with :
# snakemake -j 99 --use-conda --cluster-config cluster.json --immediate-submit --notemp --cluster 'python3 slurmSubmit.py {dependencies}'

# """
# import os
# import sys

# from snakemake.utils import read_job_properties

# jobscript = sys.argv[-1]
# job_properties = read_job_properties(jobscript)

# cmdline = "sbatch "
# for param, val in job_properties['cluster'].items():
#     cmdline += "--{param} {val} ".format(param=param, val=val)

# # Set up dependencies
# dependencies = set(sys.argv[1:-1])
# if dependencies:
#     cmdline += " --dependency=afterok:{} ".format(":".join(dependencies))

# # Adding the actual job
# cmdline += jobscript

# # remove the leading and trailing white space for the submitted jobid
# cmdline += r" | awk '{print substr($NF, 0, length($NF))}'"

# sys.stdout.write(cmdline)

# os.system(cmdline)

#!/bin/bash

#SBATCH --partition=lattice,parallel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem=1G
#SBATCH --error=metqc_sbatch_run.%J.err
#SBATCH --output=metqc_sbatch_run.%J.out

log_dir="$(pwd)"
log_file="logs/metqc-analysis.log.txt"
num_jobs=10

snakemake --cluster-config cluster.json --cluster 'sbatch --partition={cluster.partition} --cpus-per-task={cluster.cpus-per-task} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --time={cluster.time} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}' --jobs $num_jobs --use-conda &> $log_dir/$log_file

echo "finished with exit code $? at: `date`"


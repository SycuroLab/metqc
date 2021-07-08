#!/bin/bash

#SBATCH --partition=synergy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=5G
#SBATCH --error=metqc_run.%J.err
#SBATCH --output=metqc_run.%J.out

log_dir="$(pwd)"
log_file="logs/metqc-analysis.log.txt"
num_jobs=10

echo "started at: `date`"

snakemake --cluster-config cluster.json --cluster 'sbatch --partition={cluster.partition} --cpus-per-task={cluster.cpus-per-task} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --time={cluster.time} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}' --jobs $num_jobs --use-conda &> $log_dir/$log_file

output_dir=$(grep "output_dir" < config.yaml | cut -d ' ' -f2 | sed 's/"//g')
list_files=$(grep "list_files" < config.yaml | cut -d ' ' -f2 | sed 's/"//g')

snakemake_file_dir="${output_dir}/snakemake_files"
mkdir -p $snakemake_file_dir

cp $list_files $snakemake_file_dir

cp Snakefile $snakemake_file_dir
cp config.yaml $snakemake_file_dir
cp cluster.json $snakemake_file_dir
cp metqc_sbatch.sh $snakemake_file_dir 

cp -rf logs $snakemake_file_dir
cp -rf utils $snakemake_file_dir


echo "finished with exit code $? at: `date`"


# metqc

Bioinformatics pipeline for performing QC on shotgun metagenomic data.

## Overview

This pipeline is written in snakemake and designed to automate and control the submission of processes to the Synergy server at the University of Calgary. Developed by Alana Schick for the lab of Dr. Laura Sycuro. 

Input: raw paired-end fastq files.

Output: filtered and cleaned fastq files, ready to use for further downstream analysis. The pipeline also produces html reports of various QC steps. 

## Installation

To use this pipeline, clone this repository into your project directory using the following command:

```
git clone https://github.com/alanaschick/metqc.git metqc
```

Note: you need to have snakemake installed in order to run this. To install snakemake using conda, run the following line:

```
conda install -c bioconda snakemake
```

See the snakemake [webpage](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) for details.

## Config file

The pipeline requires a config file, written in yaml, to run. See the provided example file. This is the only file that should be modified before running the pipeline. Enter any custom parameters in `config.yaml`.

## Raw data and list of files.

Specify the location of your raw data files in the config file. You also need to have a list of sample names which contains the names of the samples to run the pipeline on, one sample per line. Sample names should include everything up to the R1/R2 (or 1/2) part of the file names of the raw fastq files. See `metqc_files/list_files.txt` as an example. Specify the name of your list in the config file.

## Running the pipeline on Synergy

Test the pipeline by running `snakemake -np`. This command prints out the commands to be run without actually running them. 

To run the pipeline on the synergy compute cluster, enter the following command from the project directory:

```
snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs 100 --use-conda
```
Note: the file `cluster.json` contains the parameters for the LSF job submission system that Synergy uses. They are by default the same for each process but can be modified in this file.

## Results and log files

Snakemake will create a directory for the results of the pipeline as well as a directory for log files. Logs of each step of the pipeline will be written to the `logs` directory.

## Pipeline summary

### Steps

1) QC raw reads using fastqc and multiqc. This step generates an html file called `multiqc_report_raw.html`.

2) Adapter trimming using cutadapt. This step is optional; to disable the pipeline from running cutadapt, set the parameter `run_cutadapt = FALSE` in the config file.

3) Dereplication and filtering of low complexity reads using PRINSEQ.

4) QC on filtered reads using fastqc and multiqc again. This step generates an html file called `multiqc_report_filtered.html`.

5) Remove human reads using BMtagger. Default reference for human reads is hg19 set inclusive of mitochondrial DNA from Naccache SN, Genome Research, 2014.

6) Map reads to human reference genome using BBMap. This step generates two fastq files for each sample: one of reads that mapped to the human genome and one reads that did not map to the human genome.



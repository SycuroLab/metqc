# metqc

Bioinformatics pipeline for performing QC on metagenomic data. This pipeline was developed by Alana Schick for the lab of Dr. Laura Sycuro at the University of Calgary.

## Overview

This is a robust, extensible, pipeline written in snakemake. 

## Installation

To use this pipeline, clone this repository into your project directory using the following command:

```
git clone https://github.com/alanaschick/metqc.git projectname
```

Note: you need to have snakemake installed in order to run this. To install snakemake using conda, run the following line:

```
conda install -c bioconda snakemake
```

See the snakemake [webpage](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) for details.

## Config file

The pipeline requires a config file, written in yaml, to run. See the provided example file. This is the only file that should be modified before running the pipeline. Enter any custom parameters in `config.yaml`.

## Raw data and list of files.

Ensure that all data files are in a folder called `projectname/data/rawdata/`. You also need to have a list of sample names called `ref_files/list_files.txt` which contains the names of the samples to run the pipeline on, one sample per line. Sample names should include everything up to the R1/R2 part of the file names of the raw fastq files.

## Running the pipeline

Test the pipeline by running `snakemake -np`. This command prints out the commands to be run without actually running them. 

To run the pipeline on the synergy compute cluster, enter the following command from the project directory:

```
snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs 500 --use-conda
```

## Results and log files

All output files will be placed in the `results` directory and logs of each step of the pipeline will be written to the `logs` directory.

## Pipeline summary

### Preprocessing

1) QC raw reads using fastqc and multiqc. This step generates an html file called `multiqc_report_raw.html`.

2) Adapter trimming using cutadapt. This step is optional; to disable the pipeline from running cutadapt, set the parameter `run_cutadapt = FALSE` in the config file.

3) Dereplication and filtering of low complexity reads using PRINSEQ.

4) QC on filtered reads using fastqc and multiqc again. This step generates an html file called `multiqc_report_filtered.html`.

5) Remove human reads using BMtagger. 

6) BBMap.

Notes:


QC step 3: human masking with BMtagger (software obtained from ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger and used with a human hg19 reference set inclusive of mitochondrial DNA from Naccache SN, Genome Research, 2014):
 
bmtagger.sh -b /export/snyder_work/shared/lsycuro_labshare/dbs/bmtaggerDB/hg19_rRNA_mito_Hsapiens_rna/hg19_rRNA _mito_Hsapiens_rna_reference.bitmask -x /export/snyder_work/shared/lsycuro_labshare/srprismDB/hg19_rRNA_mito_Hsapiens_rna/hg19_rRNA_mito_Hsapiens_rna_reference.srprism -T <temp directory>  -q1 -1 <forward_reads.fastq> -2 <reverse_recads.fastq)> -o  <output directory> --extract
 
QC step 4: Map reads to human reference genome and remove reads showing strong paired alignment (this gets the little bit of residual that makes it through). I currently can’t find the unadulterated reference sequence, but I’ll keep looking! Once found, we would want to make a BBMap index file. I like BBMap as a short read aligner. It will spit out some useful data files and a new readset comprised of those that do or do not meet the mapping threshold.



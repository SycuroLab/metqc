# metqc

Bioinformatics pipeline for performing QC on metagenomic data.

## Overview

This is a robust, extensible, pipeline written in snakemake. 

## Installation

To use this pipeline, clone this repository into your project directory using the following command:

```
git clone https://github.com/alanaschick/metasta.git projectname
```

Note: you need to have snakemake installed in order to run this. To install snakemake using conda, run the following line:

```
conda install -c bioconda snakemake
```

See the snakemake [webpage] (https://bitbucket.org/johanneskoester/snakemake/wiki/Home) for details.

## Config file

The pipeline requires a config file, written in yaml, to run. See the provided example file. Most options are self-explanatory and simple to setup. Store these custom parameters in `config.yaml`.

## Raw data and list of files.

Ensure that all data files are in a folder called `projectname/data/rawdata/`. 

## Running the pipeline

In order to run the pipeline, need to have a list of sample names called `list_files.txt` which contains the names of the samples to include in the analysis, one sample per line. Sample names should include everything up to the R1/R2 part of the file names of the fastq files.

Test the pipeline by running `snakemake -np`. This command prints out the commands to be run without actually running them. 

## Results

All output files will be placed in the `results` directory and logs of each command will be written to the `logs` directory.

## Pipeline summary

### Preprocessing

1) QC raw reads using fastqc and multiqc. This step generates an html file called `multiqc_report_raw.html`.

2) Possibly another QC examination step using PRINSEQ. I can’t remember the details since it’s been a while but I believe PRINSEQ generated a lot of nice optional graphics, but in the form of flat files that you had to upload to their website to view? I remember thinking it was less useful for larger sample sets, unless we wrote some code to aggregate the output files. I’m not sure if things have improved with the package? If you want to look into it a bit more and/or suggest other tools for multi-QC. We are most interested in assessing read quality across position, presence of low complexity reads, presence of residual adaptors/barcodes, presence of duplicates, and possibly presence of contaminating sequences (detection of barcode hopping and mixed human/microbial datasets). The easier it is for us to assess and visualize these attributes across 10’s-100’s of samples coming off a run, the better.

3) QC step 1: adapter/barcode trimming. This would be a step we may not always want to employ and may sometimes need to run in multiple iterations. The tool and command I have used most often is:
cutadapt (-a <adapter seq> -A <adapter seq> -m 60; Martin M, EMBnet.Journal, 2011)

4) QC step 2: deduplication and low complexity filtering with PRINSEQ (-lc_method dust -lc_threshold 7 -derep 1; Schmieder R and Edwards R Bioinformatics 2011).

5) Repeat quality evaluation with multi-FASTQC/PRINSEQ. Most important piece here is to ensure we have removed low complexity and residual end contamination/poor quality bases before human decontamination.

6) QC step 3: human masking with BMtagger (software obtained from ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger and used with a human hg19 reference set inclusive of mitochondrial DNA from Naccache SN, Genome Research, 2014):
 
bmtagger.sh -b /export/snyder_work/shared/lsycuro_labshare/dbs/bmtaggerDB/hg19_rRNA_mito_Hsapiens_rna/hg19_rRNA _mito_Hsapiens_rna_reference.bitmask -x /export/snyder_work/shared/lsycuro_labshare/srprismDB/hg19_rRNA_mito_Hsapiens_rna/hg19_rRNA_mito_Hsapiens_rna_reference.srprism -T <temp directory>  -q1 -1 <forward_reads.fastq> -2 <reverse_recads.fastq)> -o  <output directory> --extract
 
7) QC step 4: Map reads to human reference genome and remove reads showing strong paired alignment (this gets the little bit of residual that makes it through). I currently can’t find the unadulterated reference sequence, but I’ll keep looking! Once found, we would want to make a BBMap index file. I like BBMap as a short read aligner. It will spit out some useful data files and a new readset comprised of those that do or do not meet the mapping threshold.

### Assembly



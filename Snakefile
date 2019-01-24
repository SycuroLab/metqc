# **********************************
# * Snakefile for metasta pipeline *
# **********************************

# **** Variables ****

configfile: "config.yaml"

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input:
        "results/multiqc_report_raw.html"

rule fastqc_raw:
    input:
        r1 = "data/rawdata/{sample}_R1_001.fastq.gz",
        r2 = "data/rawdata/{sample}_R2_001.fastq.gz"
    output:
        r1 = "data/rawdata/fastqc/{sample}_R1_001_fastqc.html",
        r2 = "data/rawdata/fastqc/{sample}_R2_001_fastqc.html"
    conda: "envs/fastqc_env.yaml"
    shell: "fastqc -o data/rawdata/fastqc {input.r1} {input.r2}"

rule multiqc_raw:
    input:
        r1 = expand("data/rawdata/fastqc/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        r2 = expand("data/rawdata/fastqc/{sample}_R2_001_fastqc.html", sample=SAMPLES),
    output: "results/multiqc_report_raw.html"
    conda: "envs/multiqc_env.yaml"
    shell: "multiqc -f data/rawdata/fastqc -o results -n multiqc_report_raw.html"

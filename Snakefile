# **********************************
# * Snakefile for metqc pipeline *
# **********************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

import subprocess
from os.path import join

# **** Rules ****

rule all:
    input:
        expand("data/bmtagger/{sample}_bmtagged_1.fastq", sample=SAMPLES),
        expand("data/bmtagger/{sample}_bmtagged_2.fastq", sample=SAMPLES),
        "results/multiqc_report_raw.html",
        "results/multiqc_report_all.html"

rule fastqc_raw:
    input:
        r1 = join(config["path"], "{sample}_1.fastq.gz"),
        r2 = join(config["path"], "{sample}_2.fastq.gz")
    output:
        r1 = "data/fastqc_raw/{sample}_1_fastqc.html",
        r2 = "data/fastqc_raw/{sample}_2_fastqc.html"
    conda: "metqc_files/envs/fastqc_env.yaml"
    shell: "fastqc -o data/fastqc_raw {input.r1} {input.r2}"

rule multiqc_raw:
    input:
        r1 = expand("data/fastqc_raw/{sample}_1_fastqc.html", sample=SAMPLES),
        r2 = expand("data/fastqc_raw/{sample}_2_fastqc.html", sample=SAMPLES)
    output: "results/multiqc_report_raw.html"
    conda: "metqc_files/envs/multiqc_env.yaml"
    shell: "multiqc -f data/fastqc_raw -o results -n multiqc_report_raw.html"

rule cutadapt:
    input:
        r1 = join(config["path"], "{sample}_1.fastq.gz"),
        r2 = join(config["path"], "{sample}_2.fastq.gz")
    output:
        r1 = "data/trimdata/{sample}_r1_trim.fastq.gz",
        r2 = "data/trimdata/{sample}_r2_trim.fastq.gz"
    conda: "metqc_files/envs/cutadapt_env.yaml"
    shell:
            "cutadapt -m 60 --max-n {config[maxn]} -a {config[fwd_adapter]} "
            "-A {config[rev_adapter]} -o {output.r1} -p {output.r2} "
            "{input.r1} {input.r2}"

rule decompress:
    input:
        r1 = "data/trimdata/{sample}_r1_trim.fastq.gz" if config["run_cutadapt"] else join(config["path"], "{sample}_1.fastq.gz"),
        r2 = "data/trimdata/{sample}_r2_trim.fastq.gz" if config["run_cutadapt"] else join(config["path"], "{sample}_2.fastq.gz")
    output:
        r1 = "data/trimdata/{sample}_r1_trim.fastq" if config["run_cutadapt"] else join(config["path"], "{sample}_1.fastq"),
        r2 = "data/trimdata/{sample}_r2_trim.fastq" if config["run_cutadapt"] else join(config["path"], "{sample}_2.fastq")
    shell:
            "gunzip -c {input.r1} > {output.r1}; gunzip -c {input.r2} > {output.r2}"

rule prinseq:
    input:
        r1 = "data/trimdata/{sample}_r1_trim.fastq" if config["run_cutadapt"] else join(config["path"], "{sample}_1.fastq"),
        r2 = "data/trimdata/{sample}_r2_trim.fastq" if config["run_cutadapt"] else join(config["path"], "{sample}_2.fastq")
    params:
        prefix = "data/filtdata/{sample}_filtered"
    output:
        r1 = "data/filtdata/{sample}_filtered_1.fastq",
        r2 = "data/filtdata/{sample}_filtered_2.fastq"
    conda: "metqc_files/envs/prinseq_env.yaml"
    shell:
            "perl metqc_files/scripts/prinseq-lite.pl -fastq {input.r1} -fastq2 {input.r2} "
            "-trim_qual_left {config[trimleft]} -trim_qual_right {config[trimright]} "
            "-out_good {params.prefix} -out_bad null -lc_method dust -lc_threshold 7 "
            "-derep 1 -trim_qual_type {config[trim_qual_type]} -trim_qual_window "
            "{config[trim_qual_window]} -trim_qual_step {config[trim_qual_step]} "
            "-trim_qual_rule {config[trim_qual_rule]}"

rule fastqc_filt:
    input:
        r1 = "data/filtdata/{sample}_filtered_1.fastq",
        r2 = "data/filtdata/{sample}_filtered_2.fastq"
    output:
        r1 = "data/filtdata/fastqc/{sample}_filtered_1_fastqc.html",
        r2 = "data/filtdata/fastqc/{sample}_filtered_2_fastqc.html"
    conda: "metqc_files/envs/fastqc_env.yaml"
    shell: "fastqc -o data/filtdata/fastqc {input.r1} {input.r2}"

rule bmtagger:
    input:
        r1 = "data/filtdata/{sample}_filtered_1.fastq",
        r2 = "data/filtdata/{sample}_filtered_2.fastq"
    output:
        r1 = "data/bmtagger/{sample}_bmtagged_1.fastq",
        r2 = "data/bmtagger/{sample}_bmtagged_2.fastq"
    params:
        n = "data/bmtagger/{sample}_bmtagged"
    conda: "metqc_files/envs/bmtagger_env.yaml"
    shell:
        "bmtagger.sh -b {config[bmfilter_ref]} -x {config[srprism_ref]} -q 1 -1 "
        "{input.r1} -2 {input.r2} -o {params.n} -X"

rule bbmap:
    input:
        r1 = "data/filtdata/{sample}_filtered_1.fastq",
        r2 = "data/filtdata/{sample}_filtered_2.fastq"
    output:
        ur1 = "data/bbmap/{sample}_bbmapped_1.fastq",
        ur2 = "data/bbmap/{sample}_bbmapped_2.fastq",
        mr1 = "data/bbmap/{sample}_bbmapped_human_reads_1.fastq",
        mr2 = "data/bbmap/{sample}_bbmapped_human_reads_2.fastq"
    params:
        i = "data/bmtagger/{sample}_filtered_#.fastq",
        u = "data/bbmap/{sample}_bbmapped_#.fastq",
        m = "data/bbmap/{sample}_bbmapped_human_reads_#.fastq",
        pre = "{sample}"
    conda: "metqc_files/envs/bbmap_env.yaml"
    shell:
        "bbmap.sh in={params.i} outu={params.u} outm={params.m} ref={config[bbmap_ref]} nodisk scafstats=results/bbmap_stats/{params.pre}_scafstats.txt ihist=results/bbmap_stats/{params.pre}_ihist.txt statsfile=results/bbmap_stats/{params.pre}_statsfile.txt"

rule multiqc_all:
    input:
        r1 = expand("data/bbmap/{sample}_bbmapped_1.fastq", sample=SAMPLES) if config["run_bbmap"] else expand("data/filtdata/{sample}_filtered_1.fastq", sample=SAMPLES),
        r2 = expand("data/bbmap/{sample}_bbmapped_2.fastq", sample=SAMPLES) if config["run_bbmap"] else expand("data/filtdata/{sample}_filtered_2.fastq", sample=SAMPLES),
        r3 = expand("data/filtdata/fastqc/{sample}_filtered_1_fastqc.html", sample=SAMPLES),
        r4 = expand("data/filtdata/fastqc/{sample}_filtered_2_fastqc.html", sample=SAMPLES),
        r5 = expand("data/bmtagger/{sample}_bmtagged_1.fastq", sample=SAMPLES),
        r6 = expand("data/bmtagger/{sample}_bmtagged_1.fastq", sample=SAMPLES)
    output: "results/multiqc_report_all.html"
    conda: "metqc_files/envs/multiqc_env.yaml"
    shell: "multiqc . -o results -n multiqc_report_all.html -x data/fastqc_raw/"

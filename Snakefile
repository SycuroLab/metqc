# ********************************
# * Snakefile for metqc pipeline *
# ********************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Define logic ****
qc = config["qc_only"]

def all_input_reads(qc):
    if config["qc_only"]:
        return expand(config["path"]+"{sample}"+config["for"]+".fastq", sample=SAMPLES)
    else:
        if config["run_bbmap"]:
            return expand("output/bbmap/{sample}_bbmapped_1.fastq", sample=SAMPLES)
        else:
            if config["run_bmtagger"]:
                return expand("output/bmtagger/{sample}_bmtagged_1.fastq", sample=SAMPLES)
            else:
                return expand("output/prinseq/{sample}_filtered_1.fastq", sample=SAMPLES)

# **** Rules ****

rule all:
    input:
        "results/multiqc_report.html",
        "results/multiqc_report.html" if config["qc_only"] else "results/multiqc_report_prinseq_filtered.html", 
        "results/multiqc_report_bmtagger_filtered.html", all_input_reads

rule fastqc:
    input:
        r1 = config["path"]+"{sample}"+config["for"]+".fastq",
        r2 = config["path"]+"{sample}"+config["rev"]+".fastq"
    output:
        r1 = "output/fastqc/{sample}"+config["for"]+"_fastqc.html",
        r2 = "output/fastqc/{sample}"+config["rev"]+"_fastqc.html"
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o output/fastqc {input.r1} {input.r2}"

rule multiqc:
    input:
        r1 = expand("output/fastqc/{sample}"+config["for"]+"_fastqc.html", sample=SAMPLES),
        r2 = expand("output/fastqc/{sample}"+config["rev"]+"_fastqc.html", sample=SAMPLES)
    output: "results/multiqc_report.html"
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -f output/fastqc -o results -n multiqc_report.html"

rule cutadapt:
    input:
        r1 = config["path"]+"{sample}"+config["for"]+".fastq",
        r2 = config["path"]+"{sample}"+config["rev"]+".fastq"
    output:
        r1 = "output/cutadapt/{sample}_r1_trimmed.fastq",
        r2 = "output/cutadapt/{sample}_r2_trimmed.fastq"
    conda: "utils/envs/cutadapt_env.yaml"
    shell:
            "cutadapt -m {config[minlength]} --max-n {config[maxn]} -a {config[fwd_adapter]} -A {config[rev_adapter]} "
            "-o {output.r1} -p {output.r2} "
            "{input.r1} {input.r2}"

rule prinseq:
    input:
        r1 = "output/cutadapt/{sample}_r1_trimmed.fastq" if config["run_cutadapt"] else config["path"]+"{sample}"+config["for"]+".fastq",
        r2 = "output/cutadapt/{sample}_r2_trimmed.fastq" if config["run_cutadapt"] else config["path"]+"{sample}"+config["rev"]+".fastq"
    params:
        prefix = "output/prinseq/{sample}_filtered"
    output:
        r1 = "output/prinseq/{sample}_filtered_1.fastq",
        r2 = "output/prinseq/{sample}_filtered_2.fastq"
    conda: "utils/envs/prinseq_env.yaml"
    shell:
            "perl utils/scripts/prinseq-lite.pl -fastq {input.r1} -fastq2 {input.r2} "
            "-trim_left {config[trimleft]} -trim_right {config[trimright]} "
            "-out_good {params.prefix} -out_bad null -lc_method {config[lc_method]} -lc_threshold {config[lc_threshold]} "
            "-derep 1 -trim_qual_type {config[trim_qual_type]} -trim_qual_window "
            "{config[trim_qual_window]} -trim_qual_step {config[trim_qual_step]} "
            "-trim_qual_rule {config[trim_qual_rule]} -trim_qual_left {config[trim_qual_left]} "
            "-trim_qual_right {config[trim_qual_right]} -min_len {config[minlength]} "
	    "-ns_max_n {config[maxn]}"

rule fastqc_prinseq_filt:
    input:
        r1 = "output/prinseq/{sample}_filtered_1.fastq",
        r2 = "output/prinseq/{sample}_filtered_2.fastq"
    output:
        r1 = "output/prinseq/fastqc/{sample}_filtered_1_fastqc.html",
        r2 = "output/prinseq/fastqc/{sample}_filtered_2_fastqc.html"
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o output/prinseq/fastqc {input.r1} {input.r2}"

rule multiqc_prinseq_filt:
    input:
        r1 = expand("output/prinseq/fastqc/{sample}_filtered_1_fastqc.html", sample=SAMPLES),
        r2 = expand("output/prinseq/fastqc/{sample}_filtered_2_fastqc.html", sample=SAMPLES)
    output: "results/multiqc_report_prinseq_filtered.html"
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -f output/prinseq/fastqc -o results -n multiqc_report_prinseq_filtered.html"

rule bmtagger:
    input:
        r1 = "output/prinseq/{sample}_filtered_1.fastq",
        r2 = "output/prinseq/{sample}_filtered_2.fastq"
    output:
        r1 = "output/bmtagger/{sample}_bmtagged_1.fastq",
        r2 = "output/bmtagger/{sample}_bmtagged_2.fastq"
    params:
        n = "output/bmtagger/{sample}_bmtagged"
    conda: "utils/envs/bmtagger_env.yaml"
    shell:
        "bmtagger.sh -b {config[bmfilter_ref]} -x {config[srprism_ref]} -q 1 -1 "
        "{input.r1} -2 {input.r2} -o {params.n} -X"

rule fastqc_bmtagger_filt:
    input:
        r1 = "output/bmtagger/{sample}_bmtagged_1.fastq",
        r2 = "output/bmtagger/{sample}_bmtagged_2.fastq"
    output:
        r1 = "output/bmtagger/fastqc/{sample}_bmtagged_1_fastqc.html",
        r2 = "output/bmtagger/fastqc/{sample}_bmtagged_2_fastqc.html"
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o output/bmtagger/fastqc {input.r1} {input.r2}"

rule multiqc_bmtagger_filt:
    input:
        r1 = expand("output/bmtagger/fastqc/{sample}_bmtagged_1_fastqc.html", sample=SAMPLES),
        r2 = expand("output/bmtagger/fastqc/{sample}_bmtagged_2_fastqc.html", sample=SAMPLES)
    output: "results/multiqc_report_bmtagger_filtered.html"
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -f output/bmtagger/fastqc -o results -n multiqc_report_bmtagger_filtered.html"
	
rule bbmap:
    input:
        r1 = "output/bmtagger/{sample}_bmtagged_1.fastq" if config["run_bmtagger"] else "output/prinseq/{sample}_filtered_1.fastq",
        r2 = "output/bmtagger/{sample}_bmtagged_2.fastq" if config["run_bmtagger"] else "output/prinseq/{sample}_filtered_2.fastq"
    output:
        ur1 = "output/bbmap/{sample}_bbmapped_1.fastq",
        ur2 = "output/bbmap/{sample}_bbmapped_2.fastq",
        mr1 = "output/bbmap/{sample}_human_reads_1.fastq",
        mr2 = "output/bbmap/{sample}_human_reads_2.fastq"
    params:
        i = "output/bmtagger/{sample}_bmtagged_#.fastq" if config["run_bmtagger"] else "output/prinseq/{sample}_filtered_#.fastq",
        u = "output/bbmap/{sample}_bbmapped_#.fastq",
        m = "output/bbmap/{sample}_human_reads_#.fastq",
        pre = "{sample}"
    conda: "utils/envs/bbmap_env.yaml"
    shell:
        "bbmap.sh in={params.i} outu={params.u} outm={params.m} ref={config[bbmap_ref]} nodisk scafstats=results/bbmap_stats/{params.pre}_scafstats.txt ihist=results/bbmap_stats/{params.pre}_ihist.txt statsfile=results/bbmap_stats/{params.pre}_statsfile.txt"


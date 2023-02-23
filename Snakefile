# ********************************
# * Snakefile for metqc pipeline *
# ********************************

# **** Variables ****
configfile: "config.yaml"

# **** Imports ****
import pandas as pd
import os

SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# The forward read number for output files.
forward_read_num = config["forward_read_suffix"].split(".",1)[0]

# The reverse read number for output files.
reverse_read_num = config["reverse_read_suffix"].split(".",1)[0]

#reverse_read_num = os.path.splitext(config["reverse_read_suffix"])[0]
print(reverse_read_num)


## Added new code above because this fails when fastq.gz is the file type
# The forward read number for output files.
#forward_read_num = os.path.splitext(config["forward_read_suffix"])[0]

# The reverse read number for output files.
#reverse_read_num = os.path.splitext(config["reverse_read_suffix"])[0]

# **** Define logic ****
qc = config["qc_only"]

def all_input_reads(qc):
    if config["qc_only"]:
        return expand(os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]), sample=SAMPLES)
    else:
        if config["run_bbmap"]:
            return expand(os.path.join(config["output_dir"],"bbmap","{sample}_bbmapped_1.fastq"), sample=SAMPLES)
        else:
            if config["run_bmtagger"]:
                return expand(os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_1.fastq"), sample=SAMPLES)
            else:
                return expand(os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"), sample=SAMPLES)

# **** Rules ****

rule all:
    input:
        os.path.join(config["output_dir"],"multiqc","multiqc_report.html"),
        os.path.join(config["output_dir"],"multiqc","multiqc_report.html") if config["qc_only"] else os.path.join(config["output_dir"],"multiqc","multiqc_report_prinseq_filtered.html"),
        os.path.join(config["output_dir"],"multiqc","multiqc_report_bmtagger_filtered.html"), all_input_reads

rule fastqc:
    input:
        r1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        r2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        r1 = os.path.join(config["output_dir"],"fastqc","{sample}"+forward_read_num+"_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"fastqc","{sample}"+reverse_read_num+"_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"fastqc")
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc:
    input:
        r1 = expand(os.path.join(config["output_dir"],"fastqc","{sample}"+forward_read_num+"_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"fastqc","{sample}"+reverse_read_num+"_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"multiqc","multiqc_report.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"fastqc/"),
	multiqc_dir = os.path.join(config["output_dir"],"multiqc/")
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report.html"

rule cutadapt:
    input:
        r1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        r2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        r1 = os.path.join(config["output_dir"],"cutadapt","{sample}_r1_trimmed.fastq"),
        r2 = os.path.join(config["output_dir"],"cutadapt","{sample}_r2_trimmed.fastq")
    conda: "utils/envs/cutadapt_env.yaml"
    shell:
            "cutadapt -m {config[minlength]} --max-n {config[maxn]} -a {config[fwd_adapter]} -A {config[rev_adapter]} "
            "-j {config[num_cpus]} -o {output.r1} -p {output.r2} "
            "{input.r1} {input.r2}"

rule prinseq:
    input:
        r1 = os.path.join(config["output_dir"],"cutadapt","{sample}_r1_trimmed.fastq") if config["run_cutadapt"] else os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        r2 = os.path.join(config["output_dir"],"cutadapt","{sample}_r2_trimmed.fastq") if config["run_cutadapt"] else os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    params:
        prefix = os.path.join(config["output_dir"],"prinseq","{sample}_filtered")
    output:
        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
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
        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_1_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_2_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"prinseq","fastqc/")
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc_prinseq_filt:
    input:
        r1 = expand(os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_1_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_2_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"multiqc","multiqc_report_prinseq_filtered.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"prinseq","fastqc/"),
        multiqc_dir = os.path.join(config["output_dir"],"multiqc/")
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_prinseq_filtered.html"

rule bmtagger:
    input:
        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_1.fastq"),
        r2 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_2.fastq")
    params:
        n = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged")
    conda: "utils/envs/bmtagger_env.yaml"
    shell:
        "bmtagger.sh -b {config[bmfilter_ref]} -x {config[srprism_ref]} -q 1 -1 "
        "{input.r1} -2 {input.r2} -o {params.n} -X"

rule fastqc_bmtagger_filt:
    input:
        r1 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_1.fastq"),
        r2 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"bmtagger","fastqc","{sample}_bmtagged_1_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"bmtagger","fastqc","{sample}_bmtagged_2_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"bmtagger","fastqc/")
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc_bmtagger_filt:
    input:
        r1 = expand(os.path.join(config["output_dir"],"bmtagger","fastqc","{sample}_bmtagged_1_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"bmtagger","fastqc","{sample}_bmtagged_2_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"multiqc","multiqc_report_bmtagger_filtered.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"bmtagger","fastqc/"),
        multiqc_dir = os.path.join(config["output_dir"],"multiqc/")
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_bmtagger_filtered.html"
	
rule bbmap:
    input:
        r1 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_1.fastq") if config["run_bmtagger"] else os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_2.fastq") if config["run_bmtagger"] else os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
    output:
        ur1 = os.path.join(config["output_dir"],"bbmap","{sample}_bbmapped_1.fastq"),
        ur2 = os.path.join(config["output_dir"],"bbmap","{sample}_bbmapped_2.fastq"),
        mr1 = os.path.join(config["output_dir"],"bbmap","{sample}_human_reads_1.fastq"),
        mr2 = os.path.join(config["output_dir"],"bbmap","{sample}_human_reads_2.fastq")
    params:
        i = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_#.fastq") if config["run_bmtagger"] else os.path.join(config["output_dir"],"prinseq","{sample}_filtered_#.fastq"),
        u = os.path.join(config["output_dir"],"bbmap","{sample}_bbmapped_#.fastq"),
        m = os.path.join(config["output_dir"],"bbmap","{sample}_human_reads_#.fastq"),
        pre = "{sample}",
        out_dir = os.path.join(config["output_dir"],"bbmap_stats")
    conda: "utils/envs/bbmap_env.yaml"
    shell:
        "bbmap.sh in={params.i} outu={params.u} outm={params.m} ref={config[bbmap_ref]} nodisk scafstats={params.out_dir}/{params.pre}_scafstats.txt ihist={params.out_dir}/{params.pre}_ihist.txt statsfile={params.out_dir}/{params.pre}_statsfile.txt"



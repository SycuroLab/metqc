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
#print(reverse_read_num)



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
        os.path.join(config["output_dir"],"multiqc","multiqc_report_raw.html"),
        config["output_dir"]+"/qc_seqkit.csv",
        os.path.join(config["output_dir"],"multiqc","multiqc_report_raw.html") if config["qc_only"] else os.path.join(config["output_dir"],"multiqc","multiqc_report_prinseq_filtered.html"),
        os.path.join(config["output_dir"],"multiqc","multiqc_report_bmtagger_filtered.html"), all_input_reads

rule fastqc_raw:
    input:
        r1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        r2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        r1 = os.path.join(config["output_dir"],"fastqc_raw","{sample}"+forward_read_num+"_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"fastqc_raw","{sample}"+reverse_read_num+"_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"fastqc_raw")
    conda: "fastqc"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc_raw:
    input:
        r1 = expand(os.path.join(config["output_dir"],"fastqc_raw","{sample}"+forward_read_num+"_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"fastqc_raw","{sample}"+reverse_read_num+"_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"multiqc","multiqc_report_raw.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"fastqc_raw/"),
	multiqc_dir = os.path.join(config["output_dir"],"multiqc/")
    conda: "multiqc"
    shell: "multiqc -c utils/multiqc_config.yaml -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_raw.html"

rule cutadapt:
    input:
        r1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        r2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        r1 = os.path.join(config["output_dir"],"cutadapt","{sample}_r1_trimmed.fastq"),
        r2 = os.path.join(config["output_dir"],"cutadapt","{sample}_r2_trimmed.fastq")
    conda: "cutadapt"
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
    conda: "prinseq"
    shell:
            "perl utils/scripts/prinseq-lite.pl -fastq {input.r1} -fastq2 {input.r2} "
            "-trim_left {config[trimleft]} -trim_right {config[trimright]} "
            "-out_good {params.prefix} -out_bad null -lc_method {config[lc_method]} -lc_threshold {config[lc_threshold]} "
            "-derep 1 -trim_qual_type {config[trim_qual_type]} -trim_qual_window "
            "{config[trim_qual_window]} -trim_qual_step {config[trim_qual_step]} "
            "-trim_qual_rule {config[trim_qual_rule]} -trim_qual_left {config[trim_qual_left]} "
            "-trim_qual_right {config[trim_qual_right]} -min_len {config[minlength]} "
			"-ns_max_n {config[maxn]}"

rule bmtagger:
    input:
        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_1.fastq"),
        r2 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_2.fastq"),
    params:
        n = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged"),
        r3 = os.path.join(config["output_dir"],"bmtagger","bmtagger_complete.txt")
    conda: "bmtagger"
    shell:
        "bmtagger.sh -b {config[bmfilter_ref]} -x {config[srprism_ref]} -q 1 -1 {input.r1} -2 {input.r2} -o {params.n} -X;"
        " touch  {params.r3}"


rule fastqc_prinseq_filt:
    input:
        r1 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"prinseq","{sample}_filtered_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_1_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_2_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"prinseq","fastqc/")
    conda: "fastqc"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc_prinseq_filt:
    input:
        r1 = expand(os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_1_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"prinseq","fastqc","{sample}_filtered_2_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"multiqc","multiqc_report_prinseq_filtered.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"prinseq","fastqc/"),
        multiqc_dir = os.path.join(config["output_dir"],"multiqc/")
    conda: "multiqc"
    shell: "multiqc -c utils/multiqc_config.yaml -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_prinseq_filtered.html"


rule fastqc_bmtagger_filt:
    input:
        r1 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_1.fastq"),
        r2 = os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"bmtagger","fastqc","{sample}_bmtagged_1_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"bmtagger","fastqc","{sample}_bmtagged_2_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"bmtagger","fastqc/")
    conda: "fastqc"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc_bmtagger_filt:
    input:
        r1 = expand(os.path.join(config["output_dir"],"bmtagger","fastqc","{sample}_bmtagged_1_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"bmtagger","fastqc","{sample}_bmtagged_2_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"multiqc","multiqc_report_bmtagger_filtered.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"bmtagger","fastqc/"),
        multiqc_dir = os.path.join(config["output_dir"],"multiqc/")
    conda: "multiqc"
    shell: "multiqc -c utils/multiqc_config.yaml -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_bmtagger_filtered.html"
	
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
    conda: "bbmap"
    shell:
        "bbmap.sh in={params.i} outu={params.u} outm={params.m} ref={config[bbmap_ref]} nodisk scafstats={params.out_dir}/{params.pre}_scafstats.txt ihist={params.out_dir}/{params.pre}_ihist.txt statsfile={params.out_dir}/{params.pre}_statsfile.txt"



rule seqkit:
     input:
        r1 = expand(os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_1.fastq"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"bmtagger","{sample}_bmtagged_2.fastq"), sample=SAMPLES),
     output:
        raw=config["output_dir"]+"/seq_kit_raw.csv",
        prinseq=config["output_dir"] +"/seq_kit_prinseq.csv",
        bmtagger=config["output_dir"]+"/seq_kit_bmtagger.csv",
        complete=config["output_dir"]+"/seq_kit_complete.csv"
     params:
        raw=config["input_dir"],
        prinseq=directory(os.path.join(config["output_dir"],"prinseq")),
        bmtagger=directory(os.path.join(config["output_dir"],"bmtagger"))
     conda: "seqkit"
     shell:
         "seqkit stats -j {config[num_cpus]} {params.prinseq}/*_[0-9].fastq -o {output.prinseq};"
         "seqkit stats -j {config[num_cpus]} {params.bmtagger}/*.fastq -o {output.bmtagger};"
         "seqkit stats -j {config[num_cpus]} {params.raw}/*.fastq.gz -o {output.raw};"
         "touch {output.complete};"


rule host_contamination:
     input:
        raw=config["output_dir"]+"/seq_kit_raw.csv",
        prinseq=config["output_dir"] +"/seq_kit_prinseq.csv",
        bmtagger=config["output_dir"]+"/seq_kit_bmtagger.csv",
        complete=config["output_dir"]+"/seq_kit_complete.csv"
     output:
         hc=config["output_dir"]+"/qc_seqkit.csv"
     run:
        import pandas as pd
        raw=pd.read_csv(input[0],sep="\s+")
        prinseq=pd.read_csv(input[1],sep="\s+")
        bmtagger=pd.read_csv(input[2],sep="\s+")

        raw['Sample']=raw.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split(".",expand=True,).loc[:,0].str.split("_[1|2]$",expand=True,).iloc[:,0]
        raw['Type']=raw.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split(".",expand=True,).loc[:,0].str.rsplit("_",1,expand=True,).iloc[:,1]

        prinseq['Sample']=prinseq.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_filtered",expand=True,).iloc[:,0]
        prinseq['Type']=prinseq.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_filtered_",expand=True,).iloc[:,1].str.split(".fastq",expand=True,).iloc[:,0]

        bmtagger['Sample']=bmtagger.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_bmtagged",expand=True,).iloc[:,0].str.split(".fastq",expand=True,).iloc[:,0]
        bmtagger['Type']=bmtagger.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_bmtagged_",expand=True,).iloc[:,1].str.split(".fastq",expand=True,).iloc[:,0]

        raw=raw.set_index(['Sample','Type'])
        prinseq=prinseq.set_index(['Sample','Type'])
        bmtagger=bmtagger.set_index(['Sample','Type'])
        merged=bmtagger.join(prinseq,lsuffix="_bmtagged", rsuffix="_prinseq").join(raw)
        merged=merged.reset_index()
        merged.loc[:,('num_seqs','num_seqs_bmtagged','num_seqs_prinseq')] = merged.loc[:,('num_seqs','num_seqs_bmtagged','num_seqs_prinseq')].applymap(lambda src: src.replace(',', '')).astype(float)
        merged.loc[:,'Host_contamination']=(merged.loc[:,'num_seqs_prinseq']-merged.loc[:,'num_seqs_bmtagged'])
        merged.loc[:,'Clean_reads']=merged.loc[:,'num_seqs_bmtagged']
        merged.loc[:,'Low_quality_reads']=merged.loc[:,'num_seqs']-merged.loc[:,'num_seqs_prinseq']
        merged.loc[:,['Sample','Type','Low_quality_reads','Clean_reads','Host_contamination','num_seqs','num_seqs_prinseq','num_seqs_bmtagged']].to_csv(output[0],index=False)




# Host Contamination calcualtor
# This script takes the three seqkit files and generates read counts of  clean reads and host contamination 
# The script can sometimes fail at merge step, when raw R1 and R2 suffix are different from bmtagger and prinseq files. 



import pandas as pd
import os
import logging

logging.basicConfig(level=logging.DEBUG)

logging.info("The raw file R1 suffix:"+ str(snakemake.params[0]))
logging.info("The raw file R2 suffix:"+ str(snakemake.params[1]))
        
raw=pd.read_csv(snakemake.input[0],sep="\s+").astype(str)
        
prinseq=pd.read_csv(snakemake.input[1],sep="\s+").astype(str)
        
bmtagger=pd.read_csv(snakemake.input[2],sep="\s+").astype(str)

raw['Sample']=raw.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split(snakemake.params[0],expand=True).loc[:,0].str.split(snakemake.params[1],expand=True).loc[:,0]

def find_str(x):
 for p in snakemake.params:
  if p in x:
   return p.split('_')[1].replace('R','')
 return x

raw['Type']=raw.iloc[:,0].apply(find_str)

prinseq['Sample']=prinseq.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_filtered",expand=True,).iloc[:,0]
       
prinseq['Type']=prinseq.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_filtered_",expand=True,).iloc[:,1].str.split(".fastq",expand=True,).iloc[:,0]

        
bmtagger['Sample']=bmtagger.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_bmtagged",expand=True,).iloc[:,0].str.split(".fastq",expand=True,).iloc[:,0]
        
bmtagger['Type']=bmtagger.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0]).str.split("_bmtagged_",expand=True,).iloc[:,1].str.split(".fastq",expand=True,).iloc[:,0]

        
raw=raw.set_index(['Sample','Type'])

prinseq=prinseq.set_index(['Sample','Type'])
        
bmtagger=bmtagger.set_index(['Sample','Type'])

logging.warning("Make sure these three indices match in raw, prinseq and bmtagger files because they are used in the merge step:")

logging.warning(f'Raw (Sample,Type): {raw.index[0]}')
logging.warning(f'Prinseq (Sample,Type): {prinseq.index[0]}')
logging.warning(f'Bmtagger (Sample,Type): {bmtagger.index[0]}')
        
merged=bmtagger.join(prinseq,lsuffix="_bmtagged", rsuffix="_prinseq").join(raw)
        
merged=merged.reset_index()


merged.loc[:,('num_seqs','num_seqs_bmtagged','num_seqs_prinseq')] = merged.loc[:,('num_seqs','num_seqs_bmtagged','num_seqs_prinseq')].applymap(lambda x: str(x).replace(',', '')).astype(float)
        
merged.loc[:,'Host_contamination']=(merged.loc[:,'num_seqs_prinseq']-merged.loc[:,'num_seqs_bmtagged'])
        
merged.loc[:,'Clean_reads']=merged.loc[:,'num_seqs_bmtagged']
        
merged.loc[:,'Low_quality_reads']=merged.loc[:,'num_seqs']-merged.loc[:,'num_seqs_prinseq']
        
merged.loc[:,['Sample','Type','Low_quality_reads','Clean_reads','Host_contamination','num_seqs','num_seqs_prinseq','num_seqs_bmtagged']].to_csv(snakemake.output[0],index=False)

#!/usr/bin/env python3

import os
import sys
import re
import pandas as pd
import numpy as np
from itertools import takewhile
import matplotlib.pyplot as plt

print(snakemake.output[0])
print(snakemake.output[1])
print(snakemake.output[2])
print(snakemake.output[3])
print(snakemake.output[4])

merged_tables = pd.DataFrame()
merged_stats= pd.DataFrame()
for f in snakemake.input:
        with open(f) as fin:
            headers = [x.strip() for x in takewhile(lambda x: x.startswith('#'), fin)]
            names = headers[-1].split('#')[1].strip().split('\t')
            index_col = 0
        stats=pd.Series(headers[:-1]).str.split("#|\t",expand=True)
        idx=stats.iloc[:,1]
        if merged_stats.empty:
            merged_stats = stats.iloc[:,2].rename(os.path.splitext(os.path.basename(f))[0]).to_frame()
        else:
            merged_stats = pd.merge(stats.iloc[:,2].rename(os.path.splitext(os.path.basename(f))[0]).to_frame(),
                                    merged_stats,
                                    how='outer' , 
                                    left_index=True, 
                                    right_index=True
                                    )

        iIn = pd.read_csv(f, 
                          sep='\t',
                          skiprows=len(headers),
                          names = names,
                        ).fillna('')
        iIn = iIn.set_index(iIn.columns[index_col])
        if merged_tables.empty:
            merged_tables = iIn.iloc[:,0].rename(os.path.splitext(os.path.basename(f))[0]).to_frame()
        else:
            merged_tables = pd.merge(iIn.iloc[:,0].rename(os.path.splitext(os.path.basename(f))[0]).to_frame(),
                                    merged_tables,
                                    how='outer', 
                                    left_index=True, 
                                    right_index=True
                                    )
merged_stats=merged_stats.set_index(idx)
merged_stats.to_csv(snakemake.output[1])
merged_tables.fillna('0').to_csv(snakemake.output[0])
#print(merged_tables)
merged_tables.sum(axis='index').to_csv(snakemake.output[2],header=None)

#print("*****merged_stats.sum")
merged_stats.astype(np.float64).mean(axis=1).to_csv(snakemake.output[4],header=None)

## Simple plots
df=merged_stats.transpose().astype(np.float64)
fig, ax = plt.subplots(3, 3, figsize=(12,7))

m=0
for i in range(3):
     for j in range(3):
        df.hist(column = df.columns[m], bins = 20, ax=ax[i,j], figsize=(20, 18))
        if m < 4:
           m+=1
        else:
           break
     if m >=4:
        break

merged_tables.sum(axis='index').plot.hist(bins=20, ax=ax[1,2])
ax[1, 2].set_title('Reads Used per Sample')

df.plot(x="Mean", y="STDev", kind="scatter", ax=ax[2,0])
df.plot(x="Median", y="STDev", kind="scatter", ax=ax[2,1])
df.plot(x="Mode", y="STDev", kind="scatter", ax=ax[2,2])

fig.savefig(snakemake.output[3])


#argp = argparse.ArgumentParser( prog = "merge_metaphlan_tables.py",
#    description = """Performs a table join on one or more metaphlan output files.""")
#argp.add_argument( "aistms",    metavar = "input.txt", nargs = "+",
#    help = "One or more tab-delimited text tables to join" )
#argp.add_argument( '-o',    metavar = "output.txt", nargs = 1,
#    help = "Name of output file in which joined tables are saved" )

#__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" )

#argp.usage = argp.format_usage()[7:]+"\n\n\tPlease make sure to supply file paths to the files to combine. If combining 3 files (Table1.txt, Table2.txt, and Table3.txt) the call should be:\n\n\t\tpython merge_metaphlan_tables.py Table1.txt Table2.txt Table3.txt > output.txt\n\n\tA wildcard to indicate all .txt files that start with Table can be used as follows:\n\n\t\tpython merge_metaphlan_tables.py Table*.txt > output.txt"


#def main( ):
#   print("in main")
#    args = argp.parse_args( )
#    if args.o is None:
#        merge(args.aistms, sys.stdout)
#    else:
#        with open(args.o[0], 'w') as fout:
#            merge(args.aistms, fout)

#if __name__ == '__main__':
#    main()

import sys
import pandas as pd

lengths = pd.read_csv(snakemake.input.bed,sep='\t',header=None,index_col=7)[4]




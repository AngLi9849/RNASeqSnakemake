import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import os
from snakemake.shell import shell
import pandas as pd

scalefactors = pd.read_csv(snakemake.input[1], sep="\t", index_col="sample_name")

scalefactor = scalefactors.loc[str(snakemake.wildcards.sample),"scale_factor"]

bedgraph = pd.read_csv(snakemake.input[0],sep="\t",index_col=0,header=None)

bedgraph[3]=bedgraph[3]*scalefactor

bedgraph.to_csv(snakemake.output[0],sep='\t',header=None)


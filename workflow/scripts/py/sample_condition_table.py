import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

conditiontable = pd.read_table(snakemake.input[0], index_col=False, usecols=["sample_name", "condition"], header=0,)[["sample_name","condition"]]

conditiontable.to_csv(snakemake.output[0], sep="\t",index=False)

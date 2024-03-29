import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

counts = [
    pd.read_table(
        f, index_col=0, usecols=[0, 1], header=TRUE
    )
    for f in zip(snakemake.input)
]

for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix.to_csv(snakemake.output[0], sep="\t")

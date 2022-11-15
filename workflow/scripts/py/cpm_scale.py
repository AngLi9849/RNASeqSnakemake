import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


counts = [
    pd.read_table(
        f, sep="\t",index_col=0, usecols=[0, 3], header=None,
    )
    for f in zip(snakemake.input)
]

for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "chr"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()

sample_count=matrix.loc[~matrix.index.str.startswith('spikein_')]
sample_count_sum=sample_count.sum(numeric_only=True).interpolate()
sample_count_scale=pd.concat([sample_count_sum,1000000/sample_count_sum],axis=1)
sample_count_scale.columns=['sample_counts','scale_factor']
sample_count_scale.to_csv(snakemake.output[0], sep="\t")

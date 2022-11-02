import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import os
from snakemake.shell import shell
import pandas as pd

shell(
    """
    bamCoverage -b {snakemake.input.bam} -o {snakemake.output[0]} -of bedgraph --binSize {snakemake.params.bin_size} -p {snakemake.threads} -v && 
    awk -v OFS='\\t' -v name="{snakemake.wildcards.sample}" -F'\\t' '$1 !~ "spikein_" {{sum+=$4}} END {{print name, sum}}' {snakemake.output[0]} > {snakemake.output[1]} && 
    awk -v OFS='\t' -v name="{snakemake.wildcards.sample}" -F'\\t' '$1 ~ "spikein_" {{sum+=$4}} END {{print name, sum}}' {snakemake.output[0]} > {snakemake.output[2]}
    """
)

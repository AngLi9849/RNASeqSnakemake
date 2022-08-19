import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

import os
from snakemake.shell import shell
import pandas as pd

shell(
    """
    bamCoverage -b {snakemake.input.fwdbam} -o {snakemake.output.bg_fwd} -of bedgraph --binSize {snakemake.params.bin_size} -p {snakemake.threads} -v ; 
    bamCoverage -b {snakemake.input.revbam} -o {snakemake.output.bg_rev} -of bedgraph --binSize {snakemake.params.bin_size} -p {snakemake.threads} -v
    """
)

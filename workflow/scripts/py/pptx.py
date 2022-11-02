import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

# Create Powerpoint template to the configured dimensions
import pptx as pptx

prs = pptx.Presentation()
prs.slide_width = pptx.util.Cm(float(snakemake.wildcards.width))
prs.slide_height = pptx.util.Cm(float(snakemake.wildcards.height))

prs.save(snakemake.output.pptx)


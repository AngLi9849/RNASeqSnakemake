log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(dplyr)
library(tools)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(stringr)
library(officer)
library(rvg)
library(scales)


parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Initialise
source(snakemake@config[["differential_plots"]][["scripts"]][["initialise"]])
plot_list <- snakemake@input[["plots"]]

# Import wildcards as text
prefix <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["prefix"]]))
tag <- toTitleCase(as.character(snakemake@wildcards[["tag"]]))

exp <- gsub("_"," ",as.character(snakemake@wildcards[["experiment"]]))

#Add Feature section heading and
doc <- read_docx(snakemake@input[["docx"]])
doc <- body_add(doc,fpar(ftext(exp, prop=title_bold)),style = "centered")

toc <- block_toc(level=3)
tof <- block_toc(level=1,seq_id="Figure")

doc <- body_add(doc, fpar(ftext("Table of Content",prop=fp_text(bold=TRUE))))
doc <- body_add(doc,toc)
#doc <- body_add(doc,run_pagebreak())
#doc <- body_add(doc,fpar(ftext("List of Figures",prop=fp_text(bold=TRUE))))
#doc <- body_add(doc,tof)
doc <- body_add(doc,run_pagebreak())

for (i in plot_list) {

doc <- body_add(doc,block_pour_docx(i))

}

print(doc,target=snakemake@output[["docx"]])


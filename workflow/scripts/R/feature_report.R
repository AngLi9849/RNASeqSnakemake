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
valid <- toTitleCase(as.character(snakemake@wildcards[["valid"]]))
feature <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["feature"]]))
descript <- as.character(snakemake@params[["descript"]])

#Add Feature section heading and
doc <- read_docx(snakemake@input[["docx"]])
doc <- body_add(doc,fpar(ftext(feature, prop=heading_1)),style = "heading 1")
doc <- body_add(doc,fpar(ftext(feature, prop=plain)),style = "Normal")
doc <- body_add(doc,run_pagebreak())

for (i in plot_list) {

doc <- body_add(doc,block_pour_docx(i))

}

print(doc,target=snakemake@output[["docx"]])


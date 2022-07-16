sink(log)
sink(log, type="message")

# Activate libraries for differential plotting
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(stringr)
library(officer)
library(rvg)
library(scales)


# Import wildcards as characters
prefix <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["prefix"]]))
tag <- toTitleCase(as.character(snakemake@wildcards[["tag"]]))
valid <- toTitleCase(as.character(snakemake@wildcards[["valid"]]))
feature <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["feature"]]))

# Import Plotting Parameters
descript <- as.character(snakemake@params[["descript"]])

ma_n <- as.numeric(snakemake@config[["differential_plots"]][["ma_gene_name_numbers"]])
volc_n <- as.numeric(snakemake@config[["differential_plots"]][["volcano_gene_name_numbers"]])

up_col <- as.character(snakemake@config[["differential_plots"]][["up_colour"]])
down_col <- as.character(snakemake@config[["differential_plots"]][["down_colour"]])
insig_col <- as.character(snakemake@config[["differential_plots"]][["insignificant_colour"]])
background <- "white"

sig_p <- as.numeric(snakemake@config[["differential_plots"]][["significant_p"]])
undetect_p <- as.numeric(snakemake@config[["differential_plots"]][["undetect_p"]]) 

plot_dpi <- as.numeric(snakemake@config[["differential_plots"]][["dpi"]])
tick_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["axis_ticks"]])
axis_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["axis_title"]])
label_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["labels"]])
title_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["title"]])
legend_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["legend"]])
font_size <- legend_size

ppt_w <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["width"]])
ppt_h <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["height"]])
ppt_top <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["margins"]][["top"]])
ppt_bottom <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["margins"]][["bottom"]])
ppt_left <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["margins"]][["left"]])
ppt_right <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["margins"]][["right"]])

plot_w <- ppt_w - ppt_left - ppt_right
plot_h <- ppt_h - ppt_top - ppt_bottom

# Import word docx template and set fonts
doc <- read_docx(snakemake@input[["docx"]])

heading_1 <- fp_text(bold=TRUE,font.size=14)
heading_2 <- fp_text(bold=TRUE, font.size=13)
heading_3 <- fp_text(bold=TRUE,font.size=12)
plain <- fp_text(font.size=font_size)
bold <- fp_text(bold=TRUE, font.size=font_size)
italic <- fp_text(italic=TRUE, font.size=font_size)
italic_bold <- fp_text(bold=TRUE, italic=TRUE, font.size=font_size)

# Add heading to this differential analysis
heading <- gsub("_"," ",paste(exp,"vs",control, feature, analysis ,sep=" "))
doc <- body_add(doc,fpar(ftext(heading,prop=heading_2)), style="heading 2")



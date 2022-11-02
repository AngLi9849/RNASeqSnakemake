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

source(snakemake@config[["differential_plots"]][["scripts"]][["initialise"]])
doc <- read_docx(snakemake@input[["docx"]])
 

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}


# Import lfc tables and their treatment conditions

lfc <- snakemake@input[["lfc"]]
titles <- snakemake@params[["titles"]]
descripts <- snakemake@params[["descripts"]]
group_title <- as.character(snakemake@params[["group_title"]])

doc <- body_add(doc,fpar(ftext(group_title, prop=title_bold)),style = "Normal")

for (n in 1:length(titles)) {
doc <- body_add(doc,fpar(ftext(titles[[i]], prop=bold)),style = "Normal")
doc <- body_add(doc,fpar(ftext(descripts[[i]], prop=plain)),style = "Normal")
}

# Set up basic parameters
biotypes <- c(snakemake@config[["biotypes"]])

genesets <- c(snakemake@params[["genesets"]])
gene_sets <- read.csv(snakemake@config[["gene_sets"]],header=T, sep='\t', check.names=FALSE)
names(gene_sets) <- unlist(lapply(names(gene_sets),trimws))
gene_sets <- data.frame(lapply(gene_sets,trimws))
rownames(gene_sets) <- gene_sets$set_name
min_mean <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_mean_reads"]])

rankings <- c("Length","GC")

genesets

gene_sets

goi <- unlist(lapply(unlist(strsplit(paste(gene_sets$genes[(gene_sets$set_name %in% genesets) & !as.logical(gene_sets$sep)],collapse=","),",")),trimws))
sep_gene_sets <- gene_sets$set_name[(gene_sets$set_name %in% genesets) & as.logical(gene_sets$sep)]

sep_gene_sets

doc <- body_add(doc,fpar(ftext("", prop=bold)),style = "Normal")
comment <- paste("log2FoldChanges are based on at least ", min_mean, " reads across all data sets.", sep="")

doc <- body_add(doc,fpar(ftext(comment, prop=plain)),style = "Normal")
doc <- body_add(doc,run_pagebreak())
# Load in lfc data frames as a list
lfc_ls <- lapply(lfc, function(x) {read.csv(x,sep='\t',header=T,check.names=F,index="id")})
lfc_ls <- lapply(lfc, function(x) {x %>% filter(baseMean>=min_mean)})
lfc_ls

# Find common features and names
min_set <- as.logical(snakemake@config[["group_analysis"]][["min_common_set"]])
gene_df <- distinct(bind_rows(lapply(lfc_ls, function(x) {x[,c("featureID","feature_name","biotype","group")]})))
names(gene_df) <- c("featureID","feature_name","biotype","group")
head(gene_df,5)
nrow(gene_df)

#if (min_set) {
min_genes <- Reduce(intersect,lapply(lfc_ls, function(x) {x$featureID}))
gene_df <- gene_df[gene_df$featureID %in% min_genes,]
#} 

gene_n <- nrows(gene_df)
gene_n

lfc_ls <- lapply(lfc_ls,function(x) {x[x$featureID %in% gene_df$featureID,]})
lfc_ls <- lapply(lfc_ls,function(x) {x %>% mutate(log2FoldChange=ifelse(is.na(log2FoldChange),0,log2FoldChange))})
lfc_ls <- lapply(lfc_ls,function(x) {x %>% mutate(heat=(((2*(2^log2FoldChange))/(1+(2^log2FoldChange)))-1))})
# Write experiment feature difference title to each dataframe
for (n in 1:length(lfc_ls)) {
  lfc_ls[[x]]$title <- titles[[x]]
}

i_group <- c(c(""),unique(gene_df$group[gene_df$biotype %in% biotypes]),sep_gene_sets)
for (i in i_group) {

if (i =="") {
  lfc_i <- lfc_ls
} else if (i %in% gene_sets$set_name) {
  i_set_genes <- unlist(lapply(unlist(strsplit(paste(gene_sets$genes[gene_sets$set_name==i],collapse=","),",")),trimws))
  lfc_i <- lapply(lfc_ls,function(x) {x %>% filter(root_name %in% i_set_genes)})
} else if (i %in% expr$group) {
  lfc_i <- lapply(lfc_ls,function(x) {x %>% filter(group==i)})
}

i <- gsub("_"," ",i)
i_heading <- toTitleCase(ifelse(i=="","All",i))
doc <- body_add(doc,fpar(ftext(i_heading, prop=heading_1)),style = "heading 1")

for (r in rankings) {

lfc_i <- lapply(lfc_i, function(x) {
  x %>% arrange(!!sym(r)) %>% mutate(rank= 1:n()) %>% ungroup
})

doc <- body_add(doc,fpar(ftext(r, prop=heading_2)),style = "heading 2")
heat_ylab <- r
source(snakemake@config[["group_analysis"]][["scripts"]][["heat"]])
doc <- body_add(doc,fpar(ftext("Heatmap", prop=heading_3)),style = "heading 3")
doc <- body_add(doc,value=heatmap,width = 6, height = 7, res= plot_dpi,style = "centered")
doc <- body_add(doc,run_pagebreak())
}

for (n in 1:length(titles)) {

lfc_rank <- lfc_i[[n]][,c("featureID","log2FoldChange")]
names(lfc_rank) <- c("featureID","log2FoldChange")
lfc_rank <- lfc_rank %>% 
lfc_i <- 


}

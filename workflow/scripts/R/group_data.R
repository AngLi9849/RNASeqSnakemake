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
titles <- unlist(lapply(snakemake@params[["titles"]],function(x) {gsub("_"," ",as.character(x))}))
descripts <- unlist(lapply(snakemake@params[["descripts"]],function(x) {gsub("_"," ",as.character(x))}))
group_title <- gsub("_"," ",as.character(snakemake@params[["group_title"]]))

doc <- body_add(doc,fpar(ftext(group_title, prop=title_bold)),style = "Normal")

group_title
titles
descripts

heat_w <- max(min(length(titles)*1.3,6),3)

for (n in 1:length(titles)) {
doc <- body_add(doc,fpar(ftext(titles[[n]], prop=bold)),style = "Normal")
doc <- body_add(doc,fpar(ftext(descripts[[n]], prop=plain)),style = "Normal")
}

# Set up basic parameters
biotypes <- c(snakemake@config[["biotypes"]])

genesets <- c(snakemake@params[["genesets"]])
gene_sets <- read.csv(snakemake@config[["gene_sets"]],header=T, sep='\t', check.names=FALSE)
names(gene_sets) <- unlist(lapply(names(gene_sets),trimws))
gene_sets <- data.frame(lapply(gene_sets,trimws))
rownames(gene_sets) <- gene_sets$set_name
min_rpk <- as.numeric(snakemake@config[["group_analysis"]][["min_rpk"]])
min_mean <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_mean_reads"]])

rankings <- c("Length","GC")

genesets

gene_sets

goi <- unlist(lapply(unlist(strsplit(paste(gene_sets$genes[(gene_sets$set_name %in% genesets) & !as.logical(gene_sets$sep)],collapse=","),",")),trimws))
sep_gene_sets <- gene_sets$set_name[(gene_sets$set_name %in% genesets) & as.logical(gene_sets$sep)]

sep_gene_sets

doc <- body_add(doc,fpar(ftext("", prop=bold)),style = "Normal")
comment <- paste("log2 Fold Changes shown are based on at least ", min_mean, " reads and at least ", min_rpk, " reads per kilobase across all data sets.", sep="")

sig_thresh_comment <- paste("Correlations with absolute R-values of at least ", as.character(r_threshold), " are considered significant.")

sig_cor_comment <- paste("Correlated genes are highlighted in ", if (sig_cor_only) ("all") else ("only significant"), " correlations",sep="")

doc <- body_add(doc,fpar(ftext(comment, prop=plain)),style = "Normal")
doc <- body_add(doc,fpar(ftext(sig_thresh_comment, prop=plain)),style = "Normal")
doc <- body_add(doc,fpar(ftext(sig_cor_comment, prop=plain)),style = "Normal")
doc <- body_add(doc,run_pagebreak())
# Load in lfc data frames as a list
lfc_ls <- lapply(lfc, function(x) {read.csv(x,sep='\t',header=T,check.names=F,row.names=1)})
lfc_ls <- lapply(lfc_ls, function(x) {x %>% filter((RPK>=min_rpk) & (baseMean >= min_mean))})
names(lfc_ls) <- titles

# Find common features and names
min_set <- as.logical(snakemake@config[["group_analysis"]][["min_common_set"]])
gene_df <- distinct(bind_rows(lapply(lfc_ls, function(x) {x[,c("featureID","feature_name","biotype","group","exon")]})))
names(gene_df) <- c("featureID","feature_name","biotype","group","exon")
head(gene_df,5)
nrow(gene_df)

#if (min_set) {
min_genes <- Reduce(intersect,lapply(lfc_ls, function(x) {x$featureID}))
gene_df <- gene_df[gene_df$featureID %in% min_genes,]
#} 

gene_n <- nrow(gene_df)
gene_n
lfc_ls <- lapply(lfc_ls,function(x) {x[x$featureID %in% gene_df$featureID,]})
lfc_ls <- lapply(lfc_ls,function(x) {x %>% mutate(log2FoldChange=ifelse(is.na(log2FoldChange),0,log2FoldChange))})
lfc_ls <- lapply(lfc_ls,function(x) {x %>% mutate(heat=(((2*(2^log2FoldChange))/(1+(2^log2FoldChange)))-1))})

head(lfc_ls[[1]],10)


for (i in unique(gene_df$biotype[!is.na(gene_df$exon)])) {
  print(paste(i, unique(gene_df$exon[!is.na(gene_df$exon) & gene_df$biotype==i]),length(unique(gene_df$exon[!is.na(gene_df$exon) & gene_df$biotype==i]))))
  if (length(unique(gene_df$exon[!is.na(gene_df$exon) & gene_df$biotype==i]))==1) {
    gene_df$group[gene_df$biotype==i] <- gsub("_", " ", i)
  }
}

lfc_ls <- lapply(lfc_ls,function(x) {x %>% mutate(group=gene_df$group[match(featureID,gene_df$featureID)])})

# Write experiment feature difference title to each dataframe
for (n in 1:length(lfc_ls)) {
  lfc_ls[[n]]$title <- titles[[n]]
}

i_group <- c(c(""),unique(gene_df$group[gene_df$biotype %in% biotypes]),sep_gene_sets)
for (i in i_group) {

#i <- ""

if (i =="") {
  lfc_i <- lfc_ls
} else if (i %in% gene_sets$set_name) {
  i_set_genes <- unlist(lapply(unlist(strsplit(paste(gene_sets$genes[gene_sets$set_name==i],collapse=","),",")),trimws))
  lfc_i <- lapply(lfc_ls,function(x) {x %>% filter(root_name %in% i_set_genes)})
} else if (i %in% gene_df$group) {
  lfc_i <- lapply(lfc_ls,function(x) {x %>% filter(group==i)})
}

if (nrow(lfc_i[[1]]) <= 2 ) {
doc <- body_add(doc,fpar(ftext("No Qualified changes", prop=plain)),style = "Normal")
doc <- body_add(doc,run_pagebreak())

next

}

i <- gsub("_"," ",i)
i_heading <- toTitleCase(ifelse(i=="","All",i))
doc <- body_add(doc,fpar(ftext(i_heading, prop=heading_2)),style = "heading 2")

for (r in rankings) {

lfc_i <- lapply(lfc_i, function(x) {
  x %>% arrange(!!sym(r)) %>% mutate(rank= 1:n()) %>% ungroup
})
gene_n <- nrow(lfc_i[[1]])

doc <- body_add(doc,fpar(ftext(r, prop=heading_3)),style = "heading 3")
heat_ylab <- r
source(snakemake@config[["group_analysis"]][["scripts"]][["heat"]])
doc <- body_add(doc,fpar(ftext("Heatmap", prop=heading_3)),style = "heading 4")
doc <- body_add(doc,value=heatmap,width = heat_w, height = 4, res= plot_dpi,style = "centered")
doc <- body_add(doc,run_pagebreak())
}

for (n in 1:length(titles)) {

#n <- 1

main <- titles[[n]]
n_heading <- paste("Against ",main,sep=" ")
doc <- body_add(doc,fpar(ftext(n_heading, prop=heading_3)),style = "heading 3")
heat_ylab <- paste(main,"log2 Fold Change")

lfc_rank <- lfc_i[[n]][,c("featureID","log2FoldChange")]
names(lfc_rank) <- c("featureID","log2FoldChange")
lfc_rank <- lfc_rank %>% arrange(log2FoldChange) %>% mutate(rank=1:n())

lfc_i <- lapply(lfc_i, function(x) {
  x %>% mutate(rank=lfc_rank$rank[match(featureID,lfc_rank$featureID)])
})

gene_n <- nrow(lfc_rank)

source(snakemake@config[["group_analysis"]][["scripts"]][["heat"]])
doc <- body_add(doc,fpar(ftext("Heatmap", prop=heading_3)),style = "heading 4")
doc <- body_add(doc,value=heatmap,width = heat_w, height = 4, res= plot_dpi,style = "centered")
doc <- body_add(doc,run_pagebreak())
doc <- body_add(doc,fpar(ftext("Correlations", prop=heading_3)),style = "heading 4")

cor_n <- 1
for (t in (1:length(titles))[titles != main]) {
#t <- titles[2]
lfc_cor <- data.frame(lfc_i[[t]][,c("log2FoldChange","featureID","feature_name",rankings)],check.names=F)
names(lfc_cor) <- c("x_lfc","featureID","feature_name",rankings)
lfc_cor <- lfc_cor %>% mutate(y_lfc=lfc_rank$log2FoldChange[match(featureID,lfc_rank$featureID)]) %>% ungroup
source(snakemake@config[["group_analysis"]][["scripts"]][["correlation"]])
}

cor_n
cor_plots <- unlist(lapply(1:(cor_n-1), function(x) {paste("correlation_",x,sep="")}))
cor_plots_chunks <- split(cor_plots,ceiling(seq_along(cor_plots)/1))

for ( c in cor_plots_chunks) {

#c <- cor_plots_chunks[[1]]
cor_plots_ls <- lapply(c,function(x) {get(x)})
cor_height <- ceiling(length(c)/2)*5
cor_row <- ceiling(length(c)/2)

correlation_chunk <- ggarrange(plotlist=cor_plots_ls, ncol=1, nrow=cor_row)
#doc <- body_add(doc,value=correlation_1,width = 6, height = cor_height, res= plot_dpi,style = "centered")

doc <- body_add(doc,value=correlation_chunk,width = 5, height = cor_height, res= plot_dpi,style = "centered")
doc <- body_add(doc,run_pagebreak())


head(lfc_cor,20)

head(lfc_cor$Label,20)

sum(!is.na(lfc_cor$Label))

}

}

}
i_group
i_set_genes

cor_lab_n
cor_xlim
cor_ylim
cor_xbrks
cor_ybrks

length(cor_plots_ls)
cor_plots
correlation_chunk
cor_height
cor_plots_chunks
print(correlation_1)

print(doc, target = snakemake@output[["docx"]])


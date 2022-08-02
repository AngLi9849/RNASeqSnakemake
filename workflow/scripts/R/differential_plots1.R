log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(DESeq2)
library(ashr)
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

# Import snakemake parameters and inputs
dir <- as.character(snakemake@params[["dir"]])
dir.create(dir)

difference <- "expression levels"
analysis <- paste("differential", difference, "analysis")

# Import feature lengths and nucleotide content
count_length <- read.table(snakemake@input[["length"]], sep='\t',header=TRUE, check.names=FALSE)
nuc <- read.csv(snakemake@input[["nuc"]],header=T,row.names = 1, sep='\t')

# Import wildcards as text
prefix <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["prefix"]]))
tag <- toTitleCase(as.character(snakemake@wildcards[["tag"]]))
valid <- toTitleCase(as.character(snakemake@wildcards[["valid"]]))
feature <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["feature"]]))

experiment <- gsub("_"," ",as.character(snakemake@wildcards[["experiment"]]))
treatment <- as.character(snakemake@params[["treat"]])
control_cond <- as.character(snakemake@params[["control"]])
treat <-  gsub("_"," ",treatment)
control <- gsub("_"," ",control_cond)
spikein <- gsub("_"," ",as.character(snakemake@wildcards[["spikein"]]))

splice <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["splice"]]))
normaliser <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["normaliser"]]))
counting <- "read count"
counted <- gsub("_"," ",paste(tolower(splice), tolower(prefix), feature, counting, sep=" "))
norm <- gsub("_"," ",paste(spikein, normaliser, "read count", sep=" "))

biotypes <- c(snakemake@config[["biotypes"]])
goi <- c(snakemake@config[["GOI"]])

genes <- read.csv(snakemake@input[["genetab"]],sep='\t',header=F,col.names=c("gene_id","gene_name","biotype","exon_count"),check.names=F)
bed <- unique(read.csv(snakemake@input[["bed"]],sep='\t',header=F,check.names=F)[,c(7,4,8,9)])
bed[5:6] <- genes[match(bed[,2],genes$gene_id),3:4]
bed[5][is.na(bed[5])] <- bed[1][is.na(bed[5])]
genes <- bed[3:6]
names(genes) <- c("gene_id","gene_name","biotype","exon_count") 

#  Import sample config, size factors and count table and conduct DESeq2 Differential Expression Analysis
rep_pair <- as.logical(snakemake@params[["paired"]])

# Import Counts table and extract sample names involved
cts <- read.table(snakemake@input[["counts"]], sep='\t', header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)
cts <- cts[ , order(names(cts))]
samples <- names(cts)
samples

# Import size factors
size_table <- read.csv(snakemake@input[["size_table"]],header=T,sep="\t",check.names=F)
size_table <- size_table[match(samples,size_table$sample_name),]
size_factors <- as.numeric(size_table$size_factor)
names(size_factors) <- size_table$sample_name

# Calculate expression levels by rpkm
rpkm <- data.frame(lapply(names(cts), function(x) {
  cts[,paste(x)]*size_table$scale_factor[size_table$sample_name==x]
}))

names(rpkm) <- names(cts)
rownames(rpkm) <- rownames(cts)
rpkm <- rpkm/count_length$Length[match(rownames(rpkm),count_length$gene)]*(10^9)

# Convert counts table to matrix
cts_names <- row.names(cts)
cts <- sapply(cts,as.numeric)
row.names(cts) <- cts_names

# Import sample table and define ColData dataframe for deseq2
sample_table <- read.table(snakemake@params[["sample_table"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table$sample_name <- paste(sample_table$condition,"_",sample_table$protocol,"_Replicate_",sample_table$replicate,sep="")
sample_table <- sample_table[match(samples,sample_table$sample_name),]
rownames(sample_table) <- sample_table$sample_name

dds_coldata <- sample_table[,c("condition","replicate")]
coldata <- data.frame(lapply(dds_coldata,function(x) { gsub("[-_]",".",x) } ))
rownames(coldata) <- rownames(dds_coldata)

names(cts) <- rownames(coldata)
coldata
head(cts,3)

condition_col <- as.character(sample_table$colour[match(c(control_cond,treatment),sample_table$condition)])
names(condition_col) <-c(control,treat)

condition_col
# Caculate mean expression levels of each condtion
c(control_cond,treatment)
sample_table[,c("sample_name","condition")]
mean_level <- data.frame(
  lapply(c(control_cond,treatment), function(x) {
    apply(rpkm[,match(sample_table$sample_name[sample_table$condition==x],names(rpkm))],1,FUN=mean)}))
names(mean_level) <- c(control,treat)
head(mean_level,3)
write.table(rpkm, file=snakemake@output[["rpkm"]], sep='\t', row.names=F,quote=F)

compare <- list(c(control,treat))

# Calculate number of features considered in each  multi/monoexonic-biotype group
cts_genes <- data.frame(cts_names,genes$biotype[match(cts_names,genes$gene_id)],genes$exon_count[match(cts_names,genes$gene_id)])
names(cts_genes) <- c("id","biotype","exon_count")
cts_genes$exon_count <- ifelse(cts_genes$exon_count > 1, "Multiexonic", "Monoexonic")
cts_genes$group <- toTitleCase(gsub("_"," ",paste(cts_genes$exon_count,cts_genes$biotype)))
head(cts_genes,5)
# Set up DESeq2 data object
if (rep_pair){
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition + replicate)
} else {
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition)
}

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Normalise according to supplied size factors
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- size_factors
dds <- estimateDispersions(dds)
resultsNames(dds)
dds <- nbinomWaldTest(dds)
resultsNames(dds)

# Save Normalised Counts to tsv file
norm_counts <- counts(dds, normalized=T)
write.table(data.frame("gene_id"=rownames(norm_counts), norm_counts), file=snakemake@output[["normcounts"]], sep='\t', row.names=F,quote=F)

# Generate log2FoldChange shrunk results table for each experiment condition
contrast <- c("condition", gsub("[-_]",".",treatment), gsub("[-_]",".",control_cond))
res <- results(dds, contrast=contrast, parallel=parallel)
res <- lfcShrink(dds, contrast=contrast, res=res, type="ashr")

expr <- data.frame(res@listData)
rownames(expr) <- res@rownames
expr <- expr[!is.na(expr$padj),]
expr[6:8] <- genes[match(rownames(expr),genes$gene_id),2:4]
expr$exon <- ifelse(expr$exon_count>1,"Multiexonic","Monoexonic")
expr$change <- ifelse(expr$log2FoldChange>=0,"Increased","Decreased")
expr$group <- toTitleCase(gsub("_"," ",paste(expr$exon,expr$biotype)))
expr$group2<-paste(expr$change,expr$group)
expr$log10P <- -log10(expr$padj)
expr$featureID <- rownames(expr)

expr$Length <- count_length$Length[match(rownames(expr),count_length$gene)]
expr[,c("GC","AT")] <- nuc[match(rownames(expr),rownames(nuc)),c("GC","AT")]
expr$rpkm <- expr$baseMean*1000000000/expr$Length

title <- gsub("_"," ",paste(experiment, toTitleCase(analysis),sep=" "))
write.table(data.frame(expr[c(1,2,5,6:9)]),file=paste(dir,"/",title," TopTable.tsv",sep=""), sep='\t',row.names=F,quote=F)

# Initialise Plotting
source(snakemake@config[["differential_plots"]][["scripts"]][["initialise"]])
expr$colour <- ifelse(expr$padj < sig_p, ifelse(expr$log2FoldChange < 0, down_col, up_col), insig_col)

# Initialise Word Document
doc <- read_docx(snakemake@input[["docx"]])
analysis_heading <- paste( "Differential", toTitleCase(difference))
doc <- body_add(doc,fpar(ftext(analysis_heading, prop=heading_2)),style = "heading 2")

# Plot figures for features in each mono/multiexonic-biotype groups
i_group <- append(c(""),unique(expr$group[expr$biotype %in% biotypes]))
#for (i in i_group) {

i <- "Multiexonic Protein Coding"

if (i =="") {
  expr_i <- expr
} else {
  expr_i <- expr[expr$group==i,]
}

# Set file name and path
feature_i <- ifelse(i=="", feature, paste(tolower(i),feature))

file_i <- gsub("_"," ",paste(experiment, feature_i, analysis ,sep=" "))

dir_i <- paste(dir,"/",ifelse(i=="","All",i),sep="")
dir.create(dir_i)

# Write heading for this analysis group
group_heading <- gsub("_"," ",toTitleCase(ifelse(i=="","Overview",i)))
doc <- body_add(doc,fpar(ftext(group_heading, prop=heading_3)),style = "heading 3")

# Start figure counting from 1
fig_num <- run_autonum(seq_id = "Figure", pre_label = "Figure ", post_label = ". ", prop=title_bold ,tnd=3, tns="-", bkm = "plot", start_at = 1)

# Shrink Infinite log10Ps to 1.1 x maximum non-Inf log10P
max_log10P <- max(expr_i$log10P[expr_i$log10P < Inf],na.rm = T)
expr_i <-
  expr_i %>% dplyr::mutate(
    log10P = ifelse(
    log10P == Inf,
    max_log10P*1.1,
    log10P
  )
)

expr_i <- expr_i %>% arrange(change,padj) %>% group_by(change) %>% mutate(p_rank=1:n()) %>% ungroup
expr_i <- expr_i %>% arrange(change,abs(log2FoldChange)) %>% group_by(change) %>% mutate(lfc_rank=n():1) %>% ungroup

lfc_max <- max(abs(expr_i$log2FoldChange[expr_i$padj < sig_p]))
expr_i$colour <- ifelse(expr_i$padj < sig_p, ifelse(expr_i$log2FoldChange < 0, down_col, up_col), insig_col)

mean_level_i <- mean_level[rownames(mean_level) %in% expr_i$featureID,]
head(mean_level_i,5)

if (i=="") {
cts_genes_i <- cts_genes
} else {
cts_genes_i <- cts_genes[cts_genes$group==i,]
}

head(cts_genes_i,5)

title_i <- gsub("_"," ",paste(experiment, feature_i, "Differential", toTitleCase(difference),sep=" "))

# Pie Chart ===================================================
source(snakemake@config[["differential_plots"]][["scripts"]][["pie_chart"]])
# Violin Plot =================================================
source(snakemake@config[["differential_plots"]][["scripts"]][["violin_plot"]])
# MA plot ====================================================
source(snakemake@config[["differential_plots"]][["scripts"]][["ma_plot"]])

ggsave(file=paste(file_i," MA Plot.pdf",sep=""), path=dir_i,plot=ma_plot,height=9,width=12,dpi=plot_dpi)
ggsave(file=paste(file_i," MA Plot.png",sep=""), path=dir_i,plot=ma_plot,height=9,width=12,dpi=plot_dpi)

# Volcano Plot ===============================================
source(snakemake@config[["differential_plots"]][["scripts"]][["volcano_plot"]])

ggsave(file=paste(file_i," MA Plot.pdf",sep=""), path=dir_i,plot=ma_plot,height=9,width=12,dpi=plot_dpi)
ggsave(file=paste(file_i," MA Plot.png",sep=""), path=dir_i,plot=ma_plot,height=9,width=12,dpi=plot_dpi)

# Metagene ===================================================
#source(snakemake@config[["differential_plots"]][["scripts"]][["metagene"]])

# Heatmap ====================================================
#source(snakemake@config[["differential_plots"]][["scripts"]][["heatmap"]])


# GC, Length and RPKM Bias ===================================
source(snakemake@config[["differential_plots"]][["scripts"]][["bias"]])

# Summary ====================================================

sum_list <- list("pie","violin","ma","volcano")
sum_plot_list <- lapply(sum_list,function(x) {get(x)})
sum_ncol <- 2
sum_nrow <- ceiling(length(sum_plot_list)/sum_ncol)

summary <- ggarrange(plotlist=sum_plot_list,ncol=sum_ncol,nrow=sum_nrow,labels="AUTO")

summary_title <- paste(title_i, "Analysis Summary.")
summary_caption <- paste("Overviews of changes in ", feature_i, " ", difference, ". ", str_to_sentence(difference), " are compared in Deseq2 based on ", counted, " normalised to ", norm, ".", sep="")
summary_captions <- lapply(paste(sum_list,"_caption",sep=""),function(x) {get(x)})
summary_captions <- paste( "(", LETTERS[1:length(summary_captions)], "). ", summary_captions, sep="")
summary_caption <- unlist(list(summary_caption,summary_captions))

plots <- list("summary","bias","ma_plot","volcano_plot")
plot_n <- 1

# Loop for each plot listed
for ( p in plots) {
plot_p <- get(p)
title_p <- get(paste(p,"_title",sep=""))
caption_p <- get(paste(p,"_caption",sep=""))
caption_p <- format_captions(caption_p)

if (plot_n > 1) {  
fig_num <- run_autonum(seq_id = "Figure", pre_label = "Figure ", post_label = ". ", prop=title_bold ,tnd=3, tns="-", bkm = "plot")
doc <- body_add(doc,run_pagebreak())
} 


title_p <- fpar(fig_num,ftext(title_p,prop=title_prop)) 

plot_n <- plot_n + 1
if (p=="summary") {
doc <- body_add(doc,value=plot_p,width = 6, height = 8.5, res= plot_dpi,style = "centered")
doc <- body_add(doc,run_pagebreak())
} else {
doc <- body_add(doc,value=plot_p,width = 6, height = 6, res= plot_dpi,style = "centered")
}

doc <- body_add(doc,title_p)
# Loop for each line in caption
for (c in caption_p) {
doc <- body_add(doc,fpar(values=c,fp_p = fp_par(padding.top=(caption_size/2))))
}

}

#}

print(doc, target = snakemake@output[["docx"]])


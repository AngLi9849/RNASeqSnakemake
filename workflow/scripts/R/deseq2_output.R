log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(DESeq2)
library(ashr)
library(dplyr)

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

exp <- as.character(snakemake@wildcards[["experiment"]])

sample_table <- read.table(snakemake@params[["sample_table"]],header=T,sep='\t')
sample_table <- sample_table[sample_table$experiment==exp,]

a <-  as.character(snakemake@wildcards[["test"]])
control <- as.character(snakemake@wildcards[["control"]])
spikein <- as.character(snakemake@wildcards[["spikein"]])

sample_table <- sample_table[(sample_table$condition==control | sample_table$condition==a),]

biotypes <- c(snakemake@config[["biotypes"]])
goi <- c(snakemake@config[["GOI"]])

#experiments <- unique(as.character(sample_table$condition[sample_table$condition != control]))
splice <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["splice"]]))
prefix <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["prefix"]]))
tag <- toTitleCase(as.character(snakemake@wildcards[["tag"]]))
valid <- toTitleCase(as.character(snakemake@wildcards[["valid"]]))
feature <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["feature"]]))
normaliser <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["normaliser"]]))

counting <- "read count"
change <- paste(feature,counting,sep=" ")

genes <- read.csv(snakemake@input[["genetab"]],sep='\t',header=F,col.names=c("gene_id","gene_name","biotype","exon_count"),check.names=F)
bed <- unique(read.csv(snakemake@input[["bed"]],sep='\t',header=F,check.names=F)[,c(7,4,8,9)])
bed[5:6] <- genes[match(bed[,2],genes$gene_id),3:4]
bed[5][is.na(bed[5])] <- bed[1][is.na(bed[5])]

genes <- bed[3:6]
names(genes) <- c("gene_id","gene_name","biotype","exon_count") 

#  Import sample config, size factors and count table and conduct DESeq2 Differential Expression Analysis
paired <- as.logical(snakemake@params[["paired"]])

sample_table <- read.table(snakemake@params[["sample_table"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table <- sample_table[sample_table$experiment==snakemake@wildcards[["experiment"]],]
sample_table$sample_name <- paste(sample_table$condition,"_Rep_",sample_table$replicate,sep="")
rownames(sample_table) <- sample_table$sample_name

cts <- read.table(snakemake@input[["counts"]], sep='\t', header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)
cts <- cts[,match(names(cts),sample_table$sample_name)]
cts <- cts[ , order(names(cts))]
cts_names <- row.names(cts)
cts <- sapply(cts,as.numeric)
row.names(cts) <- cts_names

# Calculate number of genes considered in each  multi/monoexonic-biotype group
cts_genes <- data.frame(cts_names,genes$biotype[match(cts_names,genes$gene_id)],genes$exon_count[match(genes$gene_id)])
names(cts_genes <- c("id","biotype","exon_count")
cts_genes$exon_count <- ifelse(cts_genes$exon_count > 1, "Multiexonic", "Monoexonic")
cts_genes$group <- toTitleCase(gsub("_"," ",paste(cts_genes$exon_count,cts_genes$biotype)))
cts_genes <- data.frame(unique(cts_genes$group),lapply(unique(cts_genes$group),function(x){sum(cts_genes$group %in% x)}))

# Import size factors
size_table <- read.csv(snakemake@input[["size_table"]],header=T,sep="\t",check.names=F)
size_table <- size_table[match(size_table$sample_name,sample_table$sample_name),] 
size_factors <- as.numeric(size_table$size_factor)
names(size_factors) <- size_table$sample_name

coldata <- sample_table[,c("condition","replicate")]
coldata <- coldata[order(row.names(coldata)), , drop=F]

# Set up DESeq2 data object
if (paired){
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition + replicate)
} else {
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition)
}

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Normalise according to supplied size factors
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- size_factors
dds <- DESeq(dds,parallel=parallel)

# Save Normalised Counts to tsv file
norm_counts = counts(dds, normalized=T)
write.table(data.frame("gene_id"=rownames(norm_counts), norm_counts), file=snakemake@output[["normcounts"]], sep='\t', row.names=F,quote=F)

# Generate log2FoldChange shrunk results table for each experiment condition
contrast <- c("condition", a, control)
res <- results(dds, contrast=contrast, parallel=parallel)
res <- lfcShrink(dds, contrast=contrast, res=res, type="ashr")

expr <- data.frame(res@listData)
rownames(expr) <- res@rownames
expr <- expr[!is.na(expr$padj),]
expr[6:8] <- genes[match(rownames(expr),genes$gene_id),2:4]
expr$exon <- ifelse(expr$exon_count>1,"Multiexonic","Monoexonic")
expr$change <- ifelse(expr$log2FoldChange>=0,"Upregulated","Downregulated")
expr$group <- toTitleCase(gsub("_"," ",paste(expr$exon,expr$biotype)))
expr$group2<-paste(expr$change,expr$group)
expr$log10P <- -log10(expr$padj)

title <- gsub("_"," ",paste(a,"vs",control,toTitleCase(change),sep=" "))
write.table(data.frame(expr[c(1,2,5,6:9)]),file=paste(dir,"/",title," TopTable.tsv",sep=""), sep='\t',row.names=F,quote=F)

# Plot figures for all compared features
source(snakemake@config[["differential_plots"]][["scripts"]][["main_script"]])

print(doc, target = snakemake@output[["docx"]])


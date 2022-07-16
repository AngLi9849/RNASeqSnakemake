log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(DEXSeq)
library(ggrepel)
library(ashr)
library(dplyr)
library(tools)

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

exp <- as.character(snakemake@params[["exp"]])
paired <- snakemake@params[["paired"]]
goi <- snakemake@params[["goi"]]
sample_table <- read.table(snakemake@input[["sample_table"]],header=T,sep='\t')
control <- as.character(snakemake@params[["control"]])
experiments <- unique(as.character(sample_table$condition[sample_table$condition != control]))
feature <- gsub("([[:lower:]])([[:upper:]])",perl=TRUE,"\\1 \\2",as.character(snakemake@wildcards[["feature"]]))
change <- paste(feature,"Usage", sep=" ")
counting <- paste(feature, "Reads", sep=" ")
prefix <- gsub("([[:lower:]])([[:upper:]])",perl=TRUE,"\\1 \\2",as.character(snakemake@wildcards[["prefix"]]))
biotypes <- c(snakemake@params[["biotypes"]])

genes <- read.csv(snakemake@input[["bed"]],sep='\t',header=F,col.names=c("chr","start","end","gene_id","score","strand","feature","feature_number","variant_number","exon_count","biotype","gene_name","feature_id"),check.names=F)

up_col <- as.character(snakemake@params[["up_col"]])
down_col <- as.character(snakemake@params[["down_col"]])
p_threshold <- snakemake@params[["p_threshold"]]

ma_n <- snakemake@params[["ma_number"]]
volc_n <- snakemake@params[["volc_number"]]

cts <- read.table(snakemake@input[["counts"]], header=TRUE, sep='\t', row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)

cts <- cts[ , order(names(cts))]
cts_names <- row.names(cts)
cts <- sapply(cts,as.numeric)
row.names(cts) <- cts_names


feature_id <- row.names(cts)
group_id <- genes$gene_name[match(row.names(cts),genes$feature_id)]
length(feature_id)
length(group_id)


coldata <- sample_table[1]
for (i in (1:nrow(coldata))) {coldata$replicate[i] <- paste("Replicate ",which(row.names(sample_table)[sample_table$condition==coldata$condition[i]]==row.names(coldata)[i]),sep="")}
coldata <- coldata[order(row.names(coldata)), , drop=F]

cts_control <- cts[,sample_table$condition==control]
coldata_control <- coldata[sample_table$condition==control,]

# Loop for each expriment condition against control
for (a in experiments) {

cts_exp <- cts[,sample_table$condition==a]

cts_exp <- cbind(cts_control,cts_exp)

coldata_exp <- rbind(coldata_control,coldata[sample_table$condition==a,])

if (paired){
  full_model <- ~sample + exon + condition:exon + replicate:exon
  reduced_model <- ~sample + exon + replicate:exon
 } else {
  full_model <- ~sample + exon + condition:exon
  reduced_model <- ~sample + exon
}

dds <- DEXSeqDataSet(countData=cts_exp,sampleData=coldata_exp, featureID=feature_id, groupID=group_id, design=full_model)

rds <- estimateSizeFactors(dds)
rds <- estimateDispersions(rds,formula=full_model)
rds <- testForDEU(rds,reducedModel = reduced_model,fullModel = full_model)
rds <- estimateExonFoldChanges(rds, fitExpToVar="condition")
res <- DEXSeqResults(rds)

# Generate log2FoldChange shrunk results table for each experiment condition
title <- toTitleCase(paste(a,"vs",control,exp,change,sep=" "))
subtitle <- toTitleCase(paste(prefix, counting,sep=" "))

expr <- data.frame(res@listData[1:10])
rownames(expr) <- feature_id
colnames(expr)[10] <- "Rawlog2FoldChange"
colnames(expr)[3] <- "baseMean"
expr$gene_name <- feature_id
expr$biotype <- genes$biotype[match(feature_id,genes$feature_id)]
expr$exon_count <- genes$exon_count[match(feature_id,genes$feature_id)]
expr$exon <- ifelse(expr$exon_count>1,"Multiexonic","Monoexonic")
expr$change <- ifelse(expr$Rawlog2FoldChange>=0,"Upregulated","Downregulated")
expr$group <- toTitleCase(gsub("_"," ",paste(expr$exon,expr$biotype)))
expr$group2<-paste(expr$change,expr$group)
expr <- expr[!is.na(expr$padj),]
expr$log10P <- -log10(expr$padj)
expr$z <- qnorm(expr$pvalue)
expr$lfcSE <- abs(expr$Rawlog2FoldChange/expr$z)
ash <- ash(betahat=expr$Rawlog2FoldChange,sebetahat=expr$lfcSE, mixcompdist = "normal",method= "shrink")
expr$log2FoldChange <- ash$result$PosteriorMean
expr <- expr %>% arrange(group2,padj) %>% group_by(group2) %>% mutate(p_rank=1:n()) %>% ungroup
expr <- expr %>% arrange(group2,abs(log2FoldChange)) %>% group_by(group2) %>% mutate(lfc_rank=n():1) %>% ungroup

expr <- expr %>% arrange(change,padj) %>% group_by(change) %>% mutate(total_p_rank=1:n()) %>% ungroup
expr <- expr %>% arrange(change,abs(log2FoldChange)) %>% group_by(change) %>% mutate(total_lfc_rank=n():1) %>% ungroup

source(snakemake@params[["plot_script"]])

}


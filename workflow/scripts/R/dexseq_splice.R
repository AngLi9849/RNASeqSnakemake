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

# Set Control and Treatment Conditions
treatment <- as.character(snakemake@params[["treat"]])
control_cond <- as.character(snakemake@params[["control"]])
treat <-  gsub("_"," ",treatment)
control <- gsub("_"," ",control_cond)

# Import feature lengths and nucleotide content
express <- read.table(snakemake@input[["express"]], sep='\t',header=TRUE, check.names=FALSE)
nuc <- read.csv(snakemake@input[["nuc"]],header=T,row.names = 1, sep='\t')


# Import feature biotype and exon count information
genes <- read.csv(snakemake@input[["genetab"]],sep='\t',header=F,col.names=c("gene_id","gene_name","biotype","exon_count","chr","start","end","strand"),check.names=F)
bed <- unique(read.csv(snakemake@input[["bed"]],sep='\t',header=F,check.names=F)[,c(7,4,8,9)])
bed[5:6] <- genes[match(bed[,2],genes$gene_id),3:4]
bed[5][is.na(bed[5])] <- bed[1][is.na(bed[5])]
genes <- bed[3:6]
names(genes) <- c("gene_id","gene_name","biotype","exon_count")

#  Import sample config, size factors and count table and conduct DESeq2 Differential Expression Analysis
rep_pair <- as.logical(snakemake@params[["paired"]])

# Import size factors
size_table <- read.csv(snakemake@input[["size_table"]],header=T,sep="\t",check.names=F)
size_table <- size_table[match(samples,size_table$sample_name),]
size_factors <- as.numeric(size_table$size_factor)
names(size_factors) <- size_table$sample_name

#Import splice and unsplice counts
splice_cts <- read.table(snakemake@input[["spliced"]], header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)
unsplice_cts <- read.table(snakemake@input[["unspliced"]], header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)

splice_cts <- splice_cts[ , order(names(splice_cts))]
samples=names(splice_cts)

splice_cts_names <- row.names(splice_cts)
splice_cts <- sapply(splice_cts,as.numeric)
row.names(splice_cts) <- splice_cts_names

unsplice_cts <- unsplice_cts[ , order(names(unsplice_cts))]
unsplice_cts_names <- row.names(unsplice_cts)
unsplice_cts <- sapply(unsplice_cts,as.numeric)
row.names(unsplice_cts) <- unsplice_cts_names

feature_id <- row.names(splice_cts)
temp <- data.frame(row.names(splice_cts))
temp$group_id <- "input"

group_id <- temp$group_id

# Import sample table
sample_table <- read.table(snakemake@config[["samples"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table$sample_name <- paste(sample_table$condition,"_",sample_table$protocol,"_Replicate_",sample_table$replicate,sep="")
sample_table <- sample_table[match(samples,sample_table$sample_name),]
rownames(sample_table) <- sample_table$sample_name

feature_id <- row.names(splice_cts)
temp <- data.frame(row.names(splice_cts))
temp$group_id <- "input"

group_id <- temp$group_id

coldata <- sample_table[,c("condition","replicate")]
coldata <- coldata[order(row.names(coldata)), , drop=F]

splice_control <- splice_cts[,sample_table$condition==control_cond]
unsplice_control <- unsplice_cts[,sample_table$condition==control_cond]
coldata_control <- coldata[sample_table$condition==control_cond,]

# Loop for each expriment condition against control
#for (a in experiments) {

splice_exp <- splice_cts[,sample_table$condition==treatment]
unsplice_exp <- unsplice_cts[,sample_table$condition==treatment]

splice_cts_exp <- cbind(splice_control,splice_exp)
unsplice_cts_exp <- cbind(unsplice_control,unsplice_exp)

coldata_exp <- rbind(coldata_control,coldata[sample_table$condition==treatment,])

coldata <- data.frame(lapply(coldata_exp,function(x) { gsub("[-_]",".",x) } ))


if (rep_pair){
  full_model <- ~sample + exon + condition:exon + replicate:exon
  reduced_model <- ~sample + exon + replicate:exon
 } else {
  full_model <- ~sample + exon + condition:exon
  reduced_model <- ~sample + exon
}

dds <- DEXSeqDataSet(countData=splice_cts_exp,sampleData=coldata, featureID=feature_id, groupID=group_id, alternativeCountData = unsplice_cts_exp, design=full_model)

rds <- estimateSizeFactors(dds)
rds <- estimateDispersions(rds,formula=full_model)
rds <- testForDEU(rds,reducedModel = reduced_model,fullModel = full_model)
rds <- estimateExonFoldChanges(rds, fitExpToVar="condition")
res <- DEXSeqResults(rds)

expr <- data.frame(res@listData[1:10])
rownames(expr) <- feature_id
colnames(expr)[10] <- "Rawlog2FoldChange"
colnames(expr)[3] <- "baseMean"
#expr <- expr[!is.na(expr$padj),]
expr[11:13] <- genes[match(rownames(expr),genes$gene_id),2:4]
expr$exon <- ifelse(expr$exon_count>1,"Multiexonic","Monoexonic")
expr$change <- ifelse(expr$Rawlog2FoldChange>=0,"Upregulated","Downregulated")
expr$group <- toTitleCase(gsub("_"," ",paste(expr$exon,expr$biotype)))
expr$group2<-paste(expr$change,expr$group)
expr$featureID <- rownames(expr)
expr$Length <- express$Length[match(rownames(expr),express$featureID)]
expr[,c("GC","AT")] <- nuc[match(rownames(expr),rownames(nuc)),c("GC","AT")]
expr$rpkm <- express$rpkm[match(rownames(expr),express$featureID)]

#expr$log10P <- -log10(expr$padj)
#expr$z <- qnorm(expr$pvalue)
#expr$lfcSE <- abs(expr$Rawlog2FoldChange/expr$z)
#ash <- ash(betahat=expr$Rawlog2FoldChange,sebetahat=expr$lfcSE, mixcompdist = "normal",method= "shrink")
#expr$log2FoldChange <- ash$result$PosteriorMean
#expr <- expr %>% arrange(group2,padj) %>% group_by(group2) %>% mutate(p_rank=1:n()) %>% ungroup
#expr <- expr %>% arrange(group2,abs(log2FoldChange)) %>% group_by(group2) %>% mutate(lfc_rank=n():1) %>% ungroup

#expr <- expr %>% arrange(change,padj) %>% group_by(change) %>% mutate(total_p_rank=1:n()) %>% ungroup
#expr <- expr %>% arrange(change,abs(log2FoldChange)) %>% group_by(change) %>% mutate(total_lfc_rank=n():1) %>% ungroup

write.table(data.frame("id"=rownames(expr),expr, check.names=FALSE),file=snakemake@output[["lfc"]],sep='\t',row.names=F,quote=F)

write.table(data.frame("id"=rownames(mean_level),mean_level, check.names=FALSE),file=snakemake@output[["levels"]],sep='\t',row.names=F,quote=F)

#}


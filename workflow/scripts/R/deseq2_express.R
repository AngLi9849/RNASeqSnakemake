log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(DESeq2)
library(ashr)
library(dplyr)
library(tools)
library(stringr)


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

treat_split <- strsplit(treatment,"")
control_split <- strsplit(control_cond,"")

cond_match_ls <- match(treat_split[[1]],control_split[[1]])
cond_match_ls <- lapply(cond_match_ls,function(x) {ifelse(is.na(x),1,x)})
nonmatch <- c(Inf)
for (i in 1:length(cond_match_ls)) {if (i != cond_match_ls[i]) {nonmatch <- c(nonmatch,i)}}
cond_match <- min(nonmatch)-1

treat_str <- substr(treatment,cond_match,length(treat_split[[1]]))
control_str <- substr(control_cond,cond_match,length(control_split[[1]]))
common_str <- ifelse(cond_match<=1,"",substr(treatment,1,cond_match))

common <- gsub("_"," ",common_str)
treat <-  gsub("_"," ",treat_str)
control <- gsub("_"," ",control_str)

# Import feature lengths and nucleotide content
count_length <- read.table(snakemake@input[["length"]], sep='\t',header=TRUE, check.names=FALSE)
nuc <- read.csv(snakemake@input[["nuc"]],header=T,row.names = 1, sep='\t')

# Import feature biotype and exon count information
genes <- read.csv(snakemake@input[["genetab"]],sep='\t',header=F,col.names=c("geneID","gene_name","biotype","exon_count","chr","start","end","strand"),check.names=F)
bed <- unique(read.csv(snakemake@input[["bed"]],sep='\t',header=F,check.names=F)[,c(4,8,9,11)])
names(bed) <- c("rootID","featureID","feature_name","baseID")

genes <- genes[genes$geneID %in% bed$rootID,]
genesc[,c("rootID","featureID","feature_name","baseID")] <- bed$featureID[match(gene$geneID,bed$rootID),c("rootID","featureID","feature_name","baseID")]

#  Import sample config, size factors and count table and conduct DESeq2 Differential Expression Analysis
rep_pair <- as.logical(snakemake@params[["paired"]])

# Import Total Read count summary
total_sum <- read.table(snakemake@input[["total_sum"]], sep='\t', header=FALSE, row.names=1, check.names=FALSE)
total_sum
names(total_sum)=c("internal","spikein")


# Import Counts table and extract sample names involved
cts <- read.table(snakemake@input[["counts"]], sep='\t', header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)
cts <- cts[ , order(names(cts))]
norm_counts <- cts
samples <- names(cts)

# Import size factors
size_table <- read.csv(snakemake@input[["size_table"]],header=T,sep="\t",check.names=F)
size_table <- size_table[match(samples,size_table$sample_name),]
size_factors <- as.numeric(size_table$size_factor)
names(size_factors) <- size_table$sample_name

# Calculate expression levels by rpkm
#norm_sum <- total_sum[,c("internal")]*size_table$scale_factor[match(rownames(total_sum),size_table$sample_name)]
norm_sum <- total_sum[,c("internal")]
norm_sum <- data.frame(internal=norm_sum)
rownames(norm_sum) <- rownames(total_sum)
norm_sum

mean_sum <- mean(norm_sum$internal)

rpkm <- data.frame(lapply(names(cts), function(x) {
  cts[,paste(x)]*size_table$scale_factor[size_table$sample_name==x]/mean_sum
}))

names(rpkm) <- names(cts)
sum <- data.frame(lapply(rpkm,sum))

rownames(rpkm) <- rownames(cts)

rpkm <- rpkm/count_length$Length[match(rownames(rpkm),count_length$gene)]*(10^9)

# Convert counts table to matrix
cts_names <- row.names(cts)
cts <- sapply(cts,as.numeric)
row.names(cts) <- cts_names

# Import sample table and define ColData dataframe for deseq2
sample_table <- read.table(snakemake@config[["samples"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table$sample_name <- paste(sample_table$condition,"_",sample_table$protocol,"_Replicate_",sample_table$replicate,sep="")
sample_table <- sample_table[match(samples,sample_table$sample_name),]
rownames(sample_table) <- sample_table$sample_name

dds_coldata <- sample_table[,c("condition","replicate")]
coldata <- data.frame(lapply(dds_coldata,function(x) { gsub("[\\+|-|_]",".",x) } ))
rownames(coldata) <- rownames(dds_coldata)

names(cts) <- rownames(coldata)

# Caculate mean expression levels of each condtion
c(control_cond,treatment)
sample_table[,c("sample_name","condition")]
mean_level <- data.frame(
  lapply(c(control_cond,treatment), function(x) {
    apply(rpkm[,match(sample_table$sample_name[sample_table$condition==x],names(rpkm))],1,FUN=mean)}))
names(mean_level) <- c(control_cond,treatment)

levels <- rpkm

 

# Calculate number of features considered in each  multi/monoexonic-biotype group
cts_genes <- data.frame(cts_names,genes$biotype[match(cts_names,genes$featureID)],genes$exon_count[match(cts_names,genes$featureID)])
names(cts_genes) <- c("id","biotype","exon_count")
cts_genes$exon_count <- ifelse(cts_genes$exon_count > 1, "Multiexonic", "Monoexonic")
cts_genes$group <- toTitleCase(gsub("_"," ",paste(cts_genes$exon_count,cts_genes$biotype)))
cts_genes$sum <- rowSums(cts)
 
# Set up DESeq2 data object
if (rep_pair){
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition + replicate)
} else {
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition)
}

#dds <- dds[rowSums(counts(dds)) > 1,]

# Normalise according to supplied size factors
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- size_factors
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
resultsNames(dds)

# Save Normalised Counts to tsv file
#norm_counts <- counts(dds, normalized=T)

norm_counts <- data.frame(lapply(names(norm_counts),function(x) {norm_counts[,paste(x)]*size_table$scale_factor[size_table$sample_name==x]}),check.names=F)

names(norm_counts) <- samples
row.names(norm_counts) <- cts_names

# Generate log2FoldChange shrunk results table for each experiment condition
contrast <- c("condition", gsub("[\\+|-|_]",".",treatment), gsub("[\\+|-|_]",".",control_cond))
res <- results(dds, contrast=contrast, parallel=parallel)
#res <- lfcShrink(dds, contrast=contrast, res=res, type="ashr") # Adaptive Shrinkage messes with padj so abandoned
expr <- data.frame(res@listData)

rownames(expr) <- res@rownames
#expr <- expr[!is.na(expr$padj),]
expr[,c("featureID","feature_name","rootID","root_name","baseID","exon_count","biotype")] <- genes[match(rownames(expr),genes$featureID),c("featureID","feature_name","geneID","gene_name","baseID","exon_count","biotype")]
expr$exon <- ifelse(expr$exon_count>1,"multiexonic","monoexonic")
expr$change <- ifelse(expr$log2FoldChange>=0,"Increased","Decreased")
expr$group <- gsub("_"," ",paste(expr$exon,expr$biotype))
expr$group2<-paste(expr$change,expr$group)
#expr$log10P <- -log10(expr$padj)

#Use base feature length if lengths of tested feature is fixed
section <- snakemake@params[["section"]]
base_bed <- read.csv(snakemake@input[["base_bed"]],header=F, sep='\t', check.names=FALSE)[,c(5,8)]
names(base_bed) <- c("Length","baseID")
base_bed <- data.frame("baseID" = unique(base_bed$baseID), "Length"=lapply(unique(base_bed$baseID) function(x) {sum(base_bed$Length[base_bed$baseID==x])}), check.names=F)
use_base_length <- (section!="body" & as.logical(snakemake@params[["main_int"]]))
head(base_bed,10)
use_base_length
if (use_base_length) {
expr$Length <- base_bed$Length[match(expr$baseID,base_bed$baseID)]
} else {
expr$Length <- express$Length[match(rownames(expr),express$featureID)]
}


expr[,c("GC","AT")] <- nuc[match(rownames(expr),rownames(nuc)),c("GC","AT")]
expr$RPKM <- expr$baseMean*1000000000/(expr$Length*mean(norm_sum$internal))

write.table(data.frame("id"=rownames(rpkm),rpkm, check.names=FALSE), file=snakemake@output[["rpkm"]], sep='\t', row.names=F,quote=F,)

write.table(data.frame("id"=rownames(norm_counts), norm_counts, check.names=FALSE), file=snakemake@output[["counts"]], sep='\t', row.names=F,quote=F)

#write.table(data.frame(cts_genes[c(1,5)], check.names=FALSE),file=snakemake@output[["count_sums"]],sep='\t',row.names=F,quote=F)

write.table(data.frame("id"=rownames(expr),expr, check.names=FALSE),file=snakemake@output[["lfc"]],sep='\t',row.names=F,quote=F)

write.table(data.frame(expr[c(1,2,5,6:9)], check.names=FALSE),file=snakemake@output[["toptable"]], sep='\t',row.names=F,quote=F)

write.table(data.frame("id"=rownames(mean_level),mean_level, check.names=FALSE),file=snakemake@output[["levels"]],sep='\t',row.names=F,quote=F)




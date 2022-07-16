log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
paired <- as.logical(snakemake@params[["paired"]])

cts <- read.table(snakemake@input[["counts"]], sep='\t', header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)
cts <- cts[ , order(names(cts))]
cts_names <- row.names(cts)
cts <- sapply(cts,as.numeric)
row.names(cts) <- cts_names
spikein <- as.character(snakemake@wildcards[["spikein"]])

sample_table <- read.table(snakemake@params[["sample_table"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table <- sample_table[sample_table$experiment==snakemake@wildcards[["experiment"]],]
sample_table$sample_name <- paste(sample_table$condition,"_Rep_",sample_table$replicate,sep="")
rownames(sample_table) <- sample_table$sample_name

coldata <- sample_table[,c("condition","replicate")]
coldata <- coldata[order(row.names(coldata)), , drop=F]

if (paired){
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition + replicate)
} else {
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition)
}

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Estimate Size Factors, if spiked-in, use spikein features as controlGenes option. 
if (spikein == "spikein") {
  spikein_gene <- read.csv(file=snakemake@input[["bed"]],sep="\t")[,c(1,8)]
  spikein_gene <- spikein_gene[grepl("spikein_",spikein_gene[,1],fixed=TRUE),2]
  spikein_gene <- rownames(dds) %in% spikein_gene
  dds <- estimateSizeFactors(dds,controlGenes=spikein_gene)
} else {
  dds <- estimateSizeFactors(dds)
}
 
# Write sample size and scale factors as table
write.table(data.frame(sizeFactors(dds),1/sizeFactors(dds)),file=snakemake@output[[1]],sep='\t',row.names = TRUE, quote = FALSE, col.names = (c("sample_name\tsize_factor","scale_factor")))

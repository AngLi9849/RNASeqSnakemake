log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]], sep='\t', header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)
cts <- cts[ , order(names(cts))]
samples <- names(cts)
cts_names <- row.names(cts)
cts <- sapply(cts,as.numeric)
row.names(cts) <- cts_names
spikein <- as.character(snakemake@wildcards[["spikein"]])

sample_table <- read.table(snakemake@config[["samples"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table$sample_name <- paste(sample_table$condition,"_",sample_table$protocol,"_Replicate_",sample_table$replicate,sep="")

samples

sample_table <- sample_table[match(samples,sample_table$sample_name),]
rownames(sample_table) <- sample_table$sample_name

sample_table


coldata <- sample_table[,c("condition","replicate")]
coldata <- coldata[order(row.names(coldata)), , drop=F]

#designs <- c(~ condition + replicate, ~ condition)
#names(designs) <- c("paired", "unpaired")


#for ( n in names(designs) ) {
# Generate Scale factors for Replicate-paired analysis

head(cts,5)

coldata

designs <- ~ condition
dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=designs)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]

features <- read.csv(file=snakemake@input[["bed"]],sep="\t")[,c(1,8)]
spikein_feat <- features[grepl("spikein_",features[,1],fixed=TRUE),2]
spikein_bool <- rownames(dds) %in% spikein_feat
internal_bool <- !(rownames(dds) %in% spikein_feat)

# Estimate Size Factors, if spiked-in, use spikein features as controlGenes option.
if (spikein == "spikein") {
  dds <- estimateSizeFactors(dds,controlGenes=spikein_bool)
} else {
  dds <- estimateSizeFactors(dds,controlGenes=internal_bool)
}

# Write sample size and scale factors as table
write.table(data.frame(sizeFactors(dds),1/sizeFactors(dds)),file=snakemake@output[[1]],sep='\t',row.names = TRUE, quote = FALSE, col.names = (c("sample_name\tsize_factor","scale_factor")))
#}

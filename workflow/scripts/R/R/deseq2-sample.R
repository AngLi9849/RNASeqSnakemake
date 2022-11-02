log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE)
cts <- cts[ , order(names(cts))]

coldata <- read.table(snakemake@input[["conditions"]], header=TRUE, row.names="sample_name", check.names=FALSE)
coldata <- coldata[order(row.names(coldata)), , drop=F]

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=as.formula(snakemake@params[["model"]]))

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

# Write sample size factors and scale factors for down-stream applications
write.table(data.frame(sizeFactors(dds),1/sizeFactors(dds)),file=snakemake@output[[3]],sep='\t',row.names = TRUE, quote = FALSE, col.names = (c("sample_name\tsize_factor","scale_factor")))
#Write dds as RDS
saveRDS(dds, file=snakemake@output[[1]])
# Write normalized counts
norm_counts = counts(dds, normalized=T)
write.table(data.frame("gene"=rownames(norm_counts), norm_counts), file=snakemake@output[[2]], sep='\t', row.names=F)

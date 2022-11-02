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

spikeindds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=as.formula(snakemake@params[["model"]]))

# remove uninformative columns
spikeindds <- spikeindds[ rowSums(counts(spikeindds)) > 1, ]
# normalization and preprocessing
spikeindds <- DESeq(spikeindds, parallel=parallel)
# Write dds object as RDS
saveRDS(spikeindds, file=snakemake@output[[1]])
# Write normalized counts
spikein_norm_counts = counts(spikeindds, normalized=T)
write.table(data.frame("gene"=rownames(spikein_norm_counts), spikein_norm_counts), file=snakemake@output[[2]], sep='\t', row.names=F)

# Generate scale factors based on sample and spikein size-factors
sampledds <- readRDS(snakemake@input[["sample_rds"]])

write.table(data.frame(sizeFactors(sampledds),sizeFactors(spikeindds),sizeFactors(sampledds)/sizeFactors(spikeindds)),file=snakemake@output[[3]],sep='\t',row.names = TRUE, quote = FALSE, col.names = (c("sample_name\tsample_size_factor","spikein_size_factor","scale_factor")))

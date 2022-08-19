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

#read sample rds
alldds <- readRDS(snakemake@input[["samplerds"]])
           
#read normaliser rds from input
normaliser <- readRDS(snakemake@input[["normaliser"]])

# Introduce Spikein size factors if present
sizeFactors(alldds) <- sizeFactors(normaliser)


# remove uninformative columns
alldds <- alldds[ rowSums(counts(alldds)) > 1, ]
# normalization and preprocessing
alldds <- DESeq(alldds, parallel=parallel)


# Write dds object as RDS
saveRDS(alldds, file=snakemake@output[[1]])
# Write normalized counts
all_norm_counts = counts(alldds, normalized=T)
write.table(data.frame("gene"=rownames(all_norm_counts), all_norm_counts), file=snakemake@output[[2]], sep='\t', row.names=F)

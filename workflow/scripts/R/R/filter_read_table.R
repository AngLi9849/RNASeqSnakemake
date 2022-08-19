log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#Filter reads per gene table according to bed file
counts <- subset(read.csv(file=snakemake@input[["all_table"]],sep='\t'),read.csv(file=snakemake@input[["all_table"]],sep='\t')[[1]] %in% read.csv(file=snakemake@input[["bed"]],sep='\t')[[4]])

#Write filtered table
write.table(counts, file=snakemake@output[[1]], sep='\t', row.names=F, col.names=F)

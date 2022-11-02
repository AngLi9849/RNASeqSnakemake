library(biomaRt)
library(tidyverse)

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because 
# of very frequent "Internal Server Error"s)

bm <- useMart("ensembl")
bm <- useDataset(str_c(snakemake@params[["species"]], "_gene_ensembl"), mart=bm)

df <- read.table(snakemake@input[["counts"]], sep='\t', header=1)

g2g <- biomaRt::getBM(
            attributes = c( "ensembl_gene_id",
                            "external_gene_name",
                            "go_id"),
            filters = "ensembl_gene_id",
            values = df$gene,
            mart = bm,
            )

annotated <- merge(df, g2g, by.x="gene", by.y="ensembl_gene_id")
annotated$gene <- ifelse(annotated$external_gene_name == '', annotated$gene, annotated$external_gene_name)
annotated$external_gene_name <- NULL
write.table(annotated, snakemake@output[["symbol"]], sep='\t', row.names=F, quote=F)



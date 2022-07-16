log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(DESeq2)
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

exp <- as.character(snakemake@wildcards[["experiment"]])
analysis <- "Differential Expression Analysis"

exp <-  as.character(snakemake@wildcards[["exp"]])
control <- as.character(snakemake@wildcards[["control"]])
spikein <- as.character(snakemake@wildcards[["spikein"]])

biotypes <- c(snakemake@config[["biotypes"]])
goi <- c(snakemake@config[["GOI"]])

splice <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["splice"]]))
normaliser <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["normaliser"]]))

counting <- "read count"
change <- paste(feature,counting,sep=" ")

genes <- read.csv(snakemake@input[["genetab"]],sep='\t',header=F,col.names=c("gene_id","gene_name","biotype","exon_count"),check.names=F)
bed <- unique(read.csv(snakemake@input[["bed"]],sep='\t',header=F,check.names=F)[,c(7,4,8,9)])
bed[5:6] <- genes[match(bed[,2],genes$gene_id),3:4]
bed[5][is.na(bed[5])] <- bed[1][is.na(bed[5])]
genes <- bed[3:6]
names(genes) <- c("gene_id","gene_name","biotype","exon_count") 

#  Import sample config, size factors and count table and conduct DESeq2 Differential Expression Analysis
paired <- as.logical(snakemake@params[["paired"]])

sample_table <- read.table(snakemake@params[["sample_table"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table <- sample_table[sample_table$experiment==snakemake@wildcards[["experiment"]],]
sample_table <- sample_table[(sample_table$condition==control | sample_table$condition==a),]
sample_table$sample_name <- paste(sample_table$condition,"_Rep_",sample_table$replicate,sep="")
rownames(sample_table) <- sample_table$sample_name

coldata <- sample_table[,c("condition","replicate")]
coldata <- coldata[order(row.names(coldata)), , drop=F]

# Import size factors
size_table <- read.csv(snakemake@input[["size_table"]],header=T,sep="\t",check.names=F)
size_table <- size_table[match(size_table$sample_name,sample_table$sample_name),]
size_factors <- as.numeric(size_table$size_factor)
names(size_factors) <- size_table$sample_name

# Import feature lengths and nucleotide content
length <- read.table(snakemake@input[["length"]], sep='\t',header=TRUE, check.names=FALSE)
nuc <- read.csv(snakemake@input[["nuc"]],header=T,row.names = 1, sep='\t')

# Import Counts Table
cts <- read.table(snakemake@input[["counts"]], sep='\t', header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)
cts <- cts[,match(names(cts),sample_table$sample_name)]

rpkm <- data.frame(lapply(names(cts), function(x) {
  cts[,paste(x)]*size_table$size_factor[size_table$sample_name==x]
}))

names(rpkm) <- names(cts)
rownames(rpkm) <- rownames(cts)
rpkm <- rpkm/length$Length[match(rownames(rpkm),length$gene)]*(10^9)

cts <- cts[ , order(names(cts))]
cts_names <- row.names(cts)
cts <- sapply(cts,as.numeric)
row.names(cts) <- cts_names

# Calculate number of features considered in each  multi/monoexonic-biotype group
cts_genes <- data.frame(cts_names,genes$biotype[match(cts_names,genes$gene_id)],genes$exon_count[match(cts_names,genes$gene_id)])
names(cts_genes) <- c("id","biotype","exon_count")
cts_genes$exon_count <- ifelse(cts_genes$exon_count > 1, "Multiexonic", "Monoexonic")
cts_genes$group <- toTitleCase(gsub("_"," ",paste(cts_genes$exon_count,cts_genes$biotype)))
cts_genes <- data.frame(unique(cts_genes$group),lapply(unique(cts_genes$group),function(x){sum(cts_genes$group %in% x)}))

# Set up DESeq2 data object
if (paired){
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition + replicate)
} else {
  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition)
}

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Normalise according to supplied size factors
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- size_factors
dds <- DESeq(dds,parallel=parallel)

# Save Normalised Counts to tsv file
norm_counts <- counts(dds, normalized=T)
write.table(data.frame("gene_id"=rownames(norm_counts), norm_counts), file=snakemake@output[["normcounts"]], sep='\t', row.names=F,quote=F)

# Generate log2FoldChange shrunk results table for each experiment condition
contrast <- c("condition", exp, control)
res <- results(dds, contrast=contrast, parallel=parallel)
res <- lfcShrink(dds, contrast=contrast, res=res, type="ashr")

expr <- data.frame(res@listData)
rownames(expr) <- res@rownames
expr <- expr[!is.na(expr$padj),]
expr[6:8] <- genes[match(rownames(expr),genes$gene_id),2:4]
expr$exon <- ifelse(expr$exon_count>1,"Multiexonic","Monoexonic")
expr$change <- ifelse(expr$log2FoldChange>=0,"Upregulated","Downregulated")
expr$group <- toTitleCase(gsub("_"," ",paste(expr$exon,expr$biotype)))
expr$group2<-paste(expr$change,expr$group)
expr$log10P <- -log10(expr$padj)

expr$Length <- length$Length[match(rownames(expr),length$gene)]
expr[,c("GC","AT")] <- nuc[match(rownames(expr),rownames(nuc)),c("GC","AT")]
expr$rpkm <- expr$baseMean*1000000000/expr$Length

title <- gsub("_"," ",paste(exp,"vs",control,toTitleCase(change),sep=" "))
write.table(data.frame(expr[c(1,2,5,6:9)]),file=paste(dir,"/",title," TopTable.tsv",sep=""), sep='\t',row.names=F,quote=F)

# Initialise Plotting
source(snakemake@config[["differential_plots"]][["scripts"]][["initialise"]])
expr$colour <- ifelse(expr$padj < sig_p, ifelse(expr$log2FoldChange < 0, down_col, up_col), insig_col)

# Plot figures for features in each mono/multiexonic-biotype groups
i_group <- append(unique(expr$group[expr$biotype %in% biotypes]), "")
for (i in i_group) {

if (i =="") {
  expr_i <- expr
} else {
  expr_i <- expr[expr$group==i,]
}

# Set file name and path
file_i <- gsub("_"," ",paste(exp, "vs", control,  i, change,sep=" "))

dir_i <- paste(dir,"/",ifelse(i=="","All",i),sep="")
dir.create(dir_i)

# Set variables for titles and legends
title <- gsub("_"," ",paste(exp,"vs",control, toTitleCase(i),toTitleCase(change),sep=" "))
capt <- gsub("_"," ",paste(splice, tolower(prefix), counting, sep=" "))
norm <- gsub("_"," ",paste(spikein, tolower(normaliser), sep=" "))

# Write heading for this analysis group
group_heading <- gsub("_"," ",toTitleCase(ifelse(i=="","Overview",i)))
doc <- body_add(doc,fpar(ftext(group_heading, prop=heading_3)),style = "heading 3")

# Start figure counting from 1
fig_num <- run_autonum(seq_id = "Figure", pre_label = "Figure ", post_label = ".", prop=bold,tnd=3, tns=".", bkm = "plot", start_at = 1)

plot_n <- 1

# Shrink Infinite log10Ps to 1.1 x maximum non-Inf log10P
max_log10P <- max(expr_i$log10P[expr_i$log10P < Inf],na.rm = T)
expr_i <-
  expr_i %>% dplyr::mutate(
    log10P = ifelse(
    log10P == Inf,
    max_log10P*1.1,
    log10P
  )
)

expr_i <- expr_i %>% arrange(change,padj) %>% group_by(change) %>% mutate(p_rank=1:n()) %>% ungroup
expr_i <- expr_i %>% arrange(change,abs(log2FoldChange)) %>% group_by(change) %>% mutate(lfc_rank=n():1) %>% ungroup

max_sig_log2FC <- max(abs(expr_i$log2FoldChange[expr_i$padj < sig_p]))
expr_i$colour <- ifelse(expr_i$padj < sig_p, ifelse(expr_i$log2FoldChange < 0, down_col, up_col), insig_col)


# Pie Chart ===================================================
sig_up <- sum(expr_i$log2FoldChange[expr_i$padj < sig_p] > 0)
sig_down <- sum(expr_i$log2FoldChange[expr_i$padj < sig_p] < 0)
insig <- sum(expr_i$pad[expr_i$padj < undetect_p] >= sig_p)
undetect <- length(cts_names) - insig - sig_up - sig_down

pie_data <- data.frame(
  c("Undetectable Change","Insignificant Change","Significant Increase","Significant Decrease"),
  c(undetect,insig,sig_up,sig_down),
  c(background,insig_col,up_col,down_col)
)

names(pie_data) <- c("Category","Numbers","Colours")
pie_label <- paste(length(cts_names),"Considered")
pie_data$Category <- paste(pie_data$Numbers,pie_data$Category)

source(snakemake@config[["differential_plots"]][["scripts"]][["pie_chart"]])
pie


# Violin Plot =================================================
violin_data <- data.frame(
  unlist(
    lapply(colnames(splice_ratio_mean), function (x) { 
      splice_ratio_mean[paste(x)] 
    })),
  unlist(
    lapply(colnames(splice_ratio_mean), function (x) { 
      replicate(nrow(splice_ratio_mean),paste(x)) 
    }))
)
names(violin_data) <- c("value","condition")

source(snakemake@config[["differential_plots"]][["scripts"]][["violin"]])
violin


# MA plot ====================================================

ma_title <- paste(title, " MA Plot.", sep = "" )
ma_title <- fpar(fig_num,ftext(ma_title,prop=plain))

source(snakemake@config[["differential_plots"]][["scripts"]][["ma_plot"]])

ggsave(file=paste(file_i," MA Plot.pdf",sep=""), path=dir_i,plot=plot_i,height=9,width=12,dpi=plot_dpi)
ggsave(file=paste(file_i," MA Plot.png",sep=""), path=dir_i,plot=plot_i,height=9,width=12,dpi=plot_dpi)


# Volcano Plot ===============================================
source(snakemake@config[["differential_plots"]][["scripts"]][["volcano_plot"]])


# Metagene ===================================================
#source(snakemake@config[["differential_plots"]][["scripts"]][["metagene"]])

# Heatmap ====================================================
#source(snakemake@config[["differential_plots"]][["scripts"]][["heatmap"]])


# GC, Length and RPKM Bias ===================================
source(snakemake@config[["differential_plots"]][["scripts"]][["bias"]])
bias

# Summary ====================================================

sum_list <- list(pie,violin,ma,volcano)
sum_ncol <- 2
sum_nrow <- ceiling(length(sum_list)/sum_ncol)

summary <- ggarrange(plotlist=sum_list,ncol=sum_ncol,nrow=sum_nrow,labels="AUTO")


caption <- paste(capt, " of ", length(expr_i$gene_name), " ", toTitleCase(i), " ", feature,"s normalised by " , norm, ". ", descript, sep="")


doc <- body_add(doc,value=plot_i,style = "centered")
doc <- body_add(doc,plot_title)
doc <- body_add(doc,fpar(ftext(caption,prop=plain)))
doc <- body_add(doc,run_pagebreak())

fig_num <- run_autonum(seq_id = "Figure", pre_label = "Figure ", post_label = ".", prop=bold,tnd=3, tns=".", bkm = "plot")


}

print(doc, target = snakemake@output[["docx"]])



print(doc, target = snakemake@output[["docx"]])


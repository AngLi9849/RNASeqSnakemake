log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(dplyr)
library(tools)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(stringr)
library(officer)
library(rvg)
library(scales)


parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Import lfc table and mean levels table
expr <- read.csv(snakemake@input[["lfc"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
mean_level <- read.csv(snakemake@input[["levels"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
cts <- read.csv(snakemake@input[["counts"]],header=T,row.names = 1, sep='\t', check.names=FALSE)

# Import snakemake parameters and inputs
# Identify analysis
diff <- as.character(snakemake@wildcards[["difference"]])
difference <- as.character(diff)
analysis <- paste("differential", difference, "analysis")

if (difference == "expression") {
  difference <- "expression levels"
  difference_unit <- "expression levels (RPKM)"
} else {
  difference <- gsub("_"," ",difference)
  difference_unit <- gsub("_"," ",difference)
}

# Import wildcards as text
prefix <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["prefix"]]))
tag <- toTitleCase(as.character(snakemake@wildcards[["tag"]]))
valid <- toTitleCase(as.character(snakemake@wildcards[["valid"]]))
feature <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["feature"]]))

experiment <- gsub("_"," ",as.character(snakemake@wildcards[["experiment"]]))
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

names(mean_level) <- c(control,treat)

spikein <- gsub("_"," ",as.character(snakemake@wildcards[["spikein"]]))

contrast <- list(c(control,treat))


splice <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["splice"]]))
normaliser <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["normaliser"]]))
counting <- ifelse(diff=="splicing_ratio","splice sites read count","read count")
counted <- gsub("_"," ",paste(tolower(splice), tolower(prefix), feature, counting, sep=" "))
norm <- gsub("_"," ",paste(spikein, normaliser, "read count", sep=" "))

biotypes <- c(snakemake@config[["biotypes"]])
goi <- c(snakemake@config[["GOI"]])

#  Import sample config, size factors and count table and conduct DESeq2 Differential Expression Analysis
rep_pair <- as.logical(snakemake@params[["paired"]])

# Import sample table and define ColData dataframe for deseq2
sample_table <- read.table(snakemake@config[["samples"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table$sample_name <- paste(sample_table$condition,"_",sample_table$protocol,"_Replicate_",sample_table$replicate,sep="")
sample_table <- sample_table[sample_table$condition %in% c(control_cond, treatment),]
sample_table <- sample_table[sample_table$protocol == as.character(snakemake@params[["protocol"]]),]
rownames(sample_table) <- sample_table$sample_name

condition_col <- as.character(sample_table$colour[match(c(control_cond,treatment),sample_table$condition)])
names(condition_col) <-c(control,treat)

condition_col
compare <- list(c(control,treat))

title <- gsub("_"," ",paste(experiment, toTitleCase(analysis),sep=" "))

# Initialise Plotting
source(snakemake@config[["differential_plots"]][["scripts"]][["initialise"]])
expr$colour <- ifelse(expr$padj < sig_p, ifelse(expr$log2FoldChange < 0, down_col, up_col), insig_col)

analysis_plots <- read.csv(snakemake@config[["analysis_summary_plots"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
analysis_plots
diff
sum_list <- as.list(colnames(analysis_plots)[analysis_plots[rownames(analysis_plots)==diff,]==TRUE])
sum_list
sum_list <- c(sum_list,"ma","volcano")


# Initialise Word Document
doc <- read_docx(snakemake@input[["docx"]])
analysis_heading <- paste( "Differential", toTitleCase(difference))
doc <- body_add(doc,fpar(ftext(analysis_heading, prop=heading_2)),style = "heading 2")

min_rpkm_pc <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_rpkm_percentile"]])
# Plot figures for features in each mono/multiexonic-biotype groups
i_group <- append(c(""),unique(expr$group[expr$biotype %in% biotypes]))
for (i in i_group) {

#i <- "Multiexonic Protein Coding"

if (i =="") {
  expr_i <- expr
} else {
  expr_i <- expr[expr$group==i,]
}

if (diff=="splicing_ratio") {
splice_sum_i <- apply(cts[(cts$state=="spliced" & cts$id %in% expr_i$featureID),1:(ncol(cts)-2)],2,sum)
unsplice_sum_i <- apply(cts[(cts$state=="unspliced" & cts$id %in% expr_i$featureID),1:(ncol(cts)-2)],2,sum)
splice_sum_i
unsplice_sum_i
sum_i <- data.frame(splice_sum_i/(splice_sum_i + 2*unsplice_sum_i),check.names=F)
} else {
sum_i <- data.frame(lapply(cts[rownames(cts) %in% expr_i$featureID,],sum))
}

if (as.logical(snakemake@config[["differential_analysis"]][["use_p_adj_min_mean"]])) {
  min_mean <- max(expr$baseMean[is.na(expr$padj)])
} else {
  min_mean <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_mean_reads"]])
}


total_i <- nrow(expr_i)

min_rpkm <- quantile(expr_i$rpkm[expr_i$baseMean > 0],min_rpkm_pc/100)

insuf_i <- sum( (expr_i$baseMean < min_mean) | (expr_i$rpkm < min_rpkm) | (is.na(expr_i$pvalue)) )

expr_i <- expr_i[(expr_i$baseMean >= min_mean & expr_i$rpkm >= min_rpkm),]

expr_i$padj <- ifelse(is.na(expr_i$padj),expr_i$pvalue,expr_i$padj)
expr_i <- expr_i[!is.na(expr_i$padj),]

if (nrow(expr_i)==0) {

doc <- body_add(doc,fpar(ftext(paste("Insufficient evidence."),prop=plain)))
doc <- body_add(doc,run_pagebreak())

next

} else {

expr_i$log10P <- -log10(expr_i$padj)

mean_level_i <- mean_level[rownames(mean_level) %in% expr_i$featureID,]

# Set file name and path
feature_i <- ifelse(i=="", feature, paste(tolower(i),feature))

file_i <- gsub("_"," ",paste(experiment, feature_i, analysis ,sep=" "))

# Write heading for this analysis group
group_heading <- gsub("_"," ",toTitleCase(ifelse(i=="","All",i)))
doc <- body_add(doc,fpar(ftext(group_heading, prop=heading_3)),style = "heading 3")

# Start figure counting from 1
fig_num <- run_autonum(seq_id = "Figure", pre_label = "Figure ", post_label = ". ", prop=title_bold ,tnd=3, tns="-", bkm = "plot", start_at = 1)

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

lfc_max <- max(c(abs(expr_i$log2FoldChange[expr_i$padj < sig_p]),0))
expr_i$colour <- ifelse(expr_i$padj < sig_p, ifelse(expr_i$log2FoldChange < 0, down_col, up_col), insig_col)

title_i <- gsub("_"," ",paste(experiment, feature_i, "Differential", toTitleCase(difference),sep=" "))

# Summary Plots =============================================
for ( i in sum_list ) {{
source(snakemake@config[["differential_plots"]][["scripts"]][[paste(i)]])
}}
# Metagene ===================================================
#if (diff=="expression") {{
#source(snakemake@config[["differential_plots"]][["scripts"]][["meta"]])
#}}
# Heatmap ====================================================
#source(snakemake@config[["differential_plots"]][["scripts"]][["heat"]])


# GC, Length and RPKM Bias ===================================
source(snakemake@config[["differential_plots"]][["scripts"]][["bias"]])

# Summary ====================================================

sum_plot_list <- lapply(sum_list,function(x) {get(x)})
sum_ncol <- 2
sum_nrow <- ceiling(length(sum_plot_list)/sum_ncol)

summary <- ggarrange(plotlist=sum_plot_list,ncol=sum_ncol,nrow=sum_nrow,labels="AUTO")

summary_title <- paste(title_i, "Analysis Summary.")
summary_caption <- paste("Overviews of changes in ", feature_i, " ", difference, ". ", str_to_sentence(difference), " are compared based on ", counted, " normalised to ", norm, ".", sep="")
summary_captions <- lapply(paste(sum_list,"_caption",sep=""),function(x) {get(x)})
summary_captions <- paste( "(", LETTERS[1:length(summary_captions)], "). ", summary_captions, sep="")
summary_caption <- unlist(list(summary_caption,summary_captions))

summary_caption

plots <- list("summary","bias","ma_plot","volcano_plot")
plot_n <- 1

# Loop for each plot listed
for ( p in plots) {
plot_p <- get(p)
title_p <- get(paste(p,"_title",sep=""))
caption_p <- get(paste(p,"_caption",sep=""))
caption_p <- format_captions(caption_p)

if (plot_n > 1) {  
fig_num <- run_autonum(seq_id = "Figure", pre_label = "Figure ", post_label = ". ", prop=title_bold ,tnd=3, tns="-", bkm = "plot")
doc <- body_add(doc,run_pagebreak())
} 


title_p <- fpar(fig_num,ftext(title_p,prop=title_prop)) 

plot_n <- plot_n + 1
if (p=="summary") {
doc <- body_add(doc,value=plot_p,width = 6, height = 8.5, res= plot_dpi,style = "centered")
doc <- body_add(doc,run_pagebreak())
} else {
doc <- body_add(doc,value=plot_p,width = 6, height = 6, res= plot_dpi,style = "centered")
}

doc <- body_add(doc,title_p)
# Loop for each line in caption
for (c in caption_p) {
doc <- body_add(doc,fpar(values=c,fp_p = fp_par(padding.top=(caption_size/2))))
}

}
}
}

print(doc, target = snakemake@output[["docx"]])


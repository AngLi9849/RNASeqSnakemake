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
expr$featureID <- rownames(expr)
expr$RPKM <- expr$rpkm

mean_level <- read.csv(snakemake@input[["levels"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
cts <- read.csv(snakemake@input[["counts"]],header=T,row.names = 1, sep='\t', check.names=FALSE)

sig_bg <- read.csv(snakemake@input[["sig_bg"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
sig <- as.numeric(snakemake@params[["sig"]])
bg <- as.numeric(snakemake@params[["bg"]])
section <- snakemake@params[["section"]]
base <- gsub("_"," ",gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@params[["base"]])))
base_feat <- gsub("_"," ",gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@params[["base_feat"]])))

mx_samples <- c(snakemake@params[["samples"]])


bed <- read.csv(snakemake@input[["bed"]],header=F, sep='\t', check.names=FALSE)[,c(8,11)]
names(bed) <- c("featureID","baseID")
expr$baseID <- bed$baseID[match(expr$featureID,bed$featureID)]
head(bed,10)

base_bed <- read.csv(snakemake@input[["base_bed"]],header=F, sep='\t', check.names=FALSE)[,c(5,8)]
names(base_bed) <- c("Length","baseID")

use_base_length <- (section!="body" & as.logical(snakemake@params[["main_int"]])) 

head(base_bed,10)
use_base_length

if (use_base_length) { 

expr$Length <- base_bed$Length[match(expr$baseID,base_bed$baseID)]

}

expr$Length[1:10]

# Import snakemake parameters and inputs
# Identify analysis
diff <- as.character(snakemake@wildcards[["difference"]])
difference <- as.character(diff)
analysis <- paste("differential", difference, "analysis")

if (difference != "splicing_ratio") {

  is_antisense <- snakemake@params[["is_antisense"]]==-1
  difference <- "expression levels"
  difference_unit <- "expression levels (RPKM)"

  mx_df <- read.csv(snakemake@input[["mx_data"]],header=T, sep='\t', check.names=FALSE)
  bef_bin <- snakemake@params[["bef_bin"]]
  main_bin <- snakemake@params[["main_bin"]]
  section <- snakemake@params[["section"]]
  plot_median <- as.logical(snakemake@config[["metagene"]][["plot_median"]])
  meta_y <- ifelse(plot_median,"Median","Mean")
  meta_bin <- plotbef_bin + plotaft_bin + main_bin
  start_name <- snakemake@params[["start_name"]]
  start_name <- ifelse(start_name=="nan","Start",start_name)
  end_name <- snakemake@params[["end_name"]]
  end_name <- ifelse(end_name=="nan","End",end_name)
  if (is_antisense) {
    plotbef_len <- snakemake@params[["plotaft_len"]]
    plotaft_len <- snakemake@params[["plotbef_len"]]
    plotbef_bin <- snakemake@params[["plotaft_bin"]]
    plotaft_bin <- snakemake@params[["plotbef_bin"]]
    len_bef_n <- snakemake@params[["len_aft"]] 
    len_aft_n <- snakemake@params[["len_bef"]]
  } else {
    plotbef_len <- snakemake@params[["plotbef_len"]]
    plotaft_len <- snakemake@params[["plotaft_len"]]
    plotbef_bin <- snakemake@params[["plotbef_bin"]]
    plotaft_bin <- snakemake@params[["plotaft_bin"]]
    len_bef_n <- snakemake@params[["len_bef"]]
    len_aft_n <- snakemake@params[["len_aft"]]
  }
  meta_bin <- plotbef_bin + plotaft_bin + main_bin

plotbef_brk <- paste(ifelse(is_antisense,"-","+"), as.character(plotbef_len), sep="")
plotaft_brk <- paste(ifelse(is_antisense,"+","-"), as.character(plotaft_len), sep="")

if (section=="body") {
meta_xbrks <- c(
 if(plotbef_len>0) (0-plotbef_bin) else NULL, 
 0,
 main_bin, 
 if(plotaft_len>0) (main_bin + plotaft_bin) else NULL
)

names(meta_xbrks) <- c(
  if(plotbef_len>0) (plotbef_brk) else NULL, 
  start_name,
  end_name, 
  if(plotaft_len>0) (plotaft_brk) else NULL
)

} else {
start_name <- ifelse(start_name=="Start",paste(base_feat,toTitleCase(section)),start_name)

bef_brk_len <- signif(len_bef_n*1.5,1)/2
bef_brk <- paste(ifelse(is_antisense,"+","-"), as.character(bef_brk_len), sep="")
bef_brk_pos <- floor(signif(bef_bin*1.5,1)/2)

aft_brk_len <- signif(len_aft_n*1.5,1)/2
aft_brk <- paste(ifelse(is_antisense,"-","+"), as.character(aft_brk_len), sep="")
aft_brk_pos <- floor(signif((main_bin-bef_bin)*1.5,1)/2)

meta_xbrks <-  c(
  if(len_bef_n>0) (0-bef_brk_pos) else NULL,
  0,
  if(len_aft_n>0) (aft_brk_pos) else NULL
)

names(meta_xbrks) <- c(
  if(len_bef_n>0)(bef_brk) else NULL,
  start_name,
  if(len_aft_n>0) (aft_brk) else NULL
)

}

heat_df <- read.csv(snakemake@input[["heat_data"]],header=T,sep='\t', check.names=FALSE)
heat_x_max <- max(heat_df$Position)
heat_x_min <- min(heat_df$Position)
heat_bin <- heat_x_max - heat_x_min
heat_xbrks <- meta_xbrks*heat_bin/meta_bin 
heat_xlim <- c(heat_x_min-0.5,heat_x_max+0.5)

heat_name <- "log2FC"
heat_colours <- lapply(strsplit(as.character(snakemake@config[["heatmap"]][["heat_colours"]]),","),trimws)[[1]]
heat_lfcbrks <- c(-3:3)
heat_lfcbrks <- unlist(lapply(heat_lfcbrks, function(x) {((2^(x+1))/((2^x)+1))-1}))
names(heat_lfcbrks) <- c(-3:3)
heat_scale_pc <- as.numeric(snakemake@config[["heatmap"]][["heat_scale_pc"]])/100


heat_ranks <- c("log2FoldChange","RPKM","Length","GC")
heat_units <- c("","","bps","%")

heat_config <- data.frame(heat_ranks,heat_units)
names(heat_config) <- c("Ranking","unit")

} else {
  difference <- gsub("_"," ",difference)
  difference_unit <- gsub("_"," ",difference)
}

head(expr,10)

head(heat_df,10)

heat_bin
heat_xbrks
heat_lfcbrks
heat_colours

# Import wildcards as text
prefix <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["prefix"]]))
tag <- toTitleCase(as.character(snakemake@wildcards[["tag"]]))
valid <- toTitleCase(as.character(snakemake@wildcards[["valid"]]))
feature <- gsub("_"," ",gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["feature"]])))
title_feature <- gsub("(?<!\\w)(.)","\\U\\1", feature, perl = TRUE)

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

genesets <- c(snakemake@params[["genesets"]])
gene_sets <- read.csv(snakemake@config[["gene_sets"]],header=T, sep='\t', check.names=FALSE)
names(gene_sets) <- unlist(lapply(names(gene_sets),trimws))
gene_sets <- data.frame(lapply(gene_sets,trimws))
rownames(gene_sets) <- gene_sets$set_name

goi <- unlist(lapply(unlist(strsplit(paste(gene_sets[genesets,]$genes,collapse=","),",")),trimws))

#goi <- c(snakemake@params[["GOI"]])

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

sample_table$cond_abbr <- c(control,treat)[match(sample_table$condition, c(control_cond,treatment))]
sample_table$sample_abbr <- paste(sample_table$cond_abbr," Rep",sample_table$replicate,sep="")

condition_col
compare <- list(c(control,treat))

title <- gsub("_"," ",paste(experiment, toTitleCase(analysis),sep=" "))

# Initialise Plotting
source(snakemake@config[["differential_plots"]][["scripts"]][["initialise"]])
expr$colour <- ifelse(expr$padj < sig_p, ifelse(expr$log2FoldChange < 0, down_col, up_col), insig_col)

head(expr,10)

for (i in unique(expr$biotype[!is.na(expr$exon)])) {
  print(paste(i, unique(expr$exon[!is.na(expr$exon) & expr$biotype==i]),length(unique(expr$exon[!is.na(expr$exon) & expr$biotype==i]))))
  if (length(unique(expr$exon[!is.na(expr$exon) & expr$biotype==i]))==1) {
    expr$group[expr$biotype==i] <- gsub("_", " ", i)
  }
}

head(expr[expr$biotype=="snRNA",],5)

analysis_plots <- read.csv(snakemake@config[["analysis_summary_plots"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
analysis_plots
sum_list <- as.list(colnames(analysis_plots)[analysis_plots[rownames(analysis_plots)==diff,]==TRUE])
sum_list
sum_list <- c(sum_list,"ma","volcano")


# Initialise Word Document
doc <- read_docx(snakemake@input[["docx"]])
min_rpkm_pc <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_rpkm_percentile"]])

# Plot figures for features in each mono/multiexonic-biotype groups
i_group <- append(c(""),unique(expr$group[expr$biotype %in% biotypes]))
for (i in i_group) {

#i <- "multiexonic protein coding"

if (i =="") {
  expr_i <- expr
} else {
  expr_i <- expr[expr$group==i,]
}

head(expr_i, 5)

if (diff=="splicing_ratio") {
splice_sum_i <- apply(cts[(cts$state=="spliced" & cts$id %in% expr_i$featureID),1:(ncol(cts)-2)],2,sum)
unsplice_sum_i <- apply(cts[(cts$state=="unspliced" & cts$id %in% expr_i$featureID),1:(ncol(cts)-2)],2,sum)
splice_sum_i
unsplice_sum_i
sum_i <- data.frame(splice_sum_i/(splice_sum_i + 2*unsplice_sum_i),check.names=F)
} else {
sum_i <- data.frame(lapply(cts[rownames(cts) %in% expr_i$featureID,],sum),check.names=F)
sum_i <- t(sum_i)
}
head(sum_i,5)

if (as.logical(snakemake@config[["differential_analysis"]][["use_p_adj_min_mean"]])) {
  min_mean <- max(expr_i$baseMean[is.na(expr_i$padj)])
} else {
  min_mean <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_mean_reads"]])
}

config_min_mean <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_mean_reads"]])

heat_min_reads <- snakemake@config[["heatmap"]][["min_reads"]]
expr_heat <- expr_i[expr_i$baseMean >= heat_min_reads,]

total_i <- nrow(expr_i)

min_rpkm <- quantile(expr_i$rpkm[expr_i$baseMean > 0],min_rpkm_pc/100)
meta_gene_n <- sum(rownames(sig_bg)[ (sig_bg$sig2bg >= sig & sig_bg$bg2sig >= bg) ] %in% expr_i$featureID[expr_i$baseMean >= config_min_mean & expr_i$baseMean >= min_rpkm])



insuf_i <- sum( (expr_i$baseMean < min_mean) | (expr_i$rpkm < min_rpkm) | (is.na(expr_i$pvalue)) )

expr_i <- expr_i[(expr_i$baseMean >= min_mean & expr_i$rpkm >= min_rpkm),]

expr_i$padj <- ifelse(is.na(expr_i$padj),expr_i$pvalue,expr_i$padj)
expr_i <- expr_i[!is.na(expr_i$padj),]


feature_i <- ifelse(i=="", feature, paste(i,feature))
base_i <- ifelse(i=="", paste("base", base), paste("base", i,base)) 

title_feature_i <- gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE)
title_base_i <- gsub("(?<!\\w)(.)","\\U\\1", base_i, perl = TRUE)


# Write heading for this analysis group
group_heading <- gsub("_"," ",toTitleCase(ifelse(i=="","All",i)))
doc <- body_add(doc,fpar(ftext(group_heading, prop=heading_3)),style = "heading 3")

group_label <- gsub("exonic ","exonic\n",group_heading)

head(expr_i,5)
if (nrow(expr_i)==0) {

doc <- body_add(doc,fpar(ftext(paste("Insufficient evidence."),prop=plain)))
doc <- body_add(doc,run_pagebreak())

next

} else {

expr_i$log10P <- -log10(expr_i$padj)

mean_level_i <- mean_level[rownames(mean_level) %in% expr_i$featureID,]

# Set file name and path
file_i <- gsub("_"," ",paste(experiment, feature_i, analysis ,sep=" "))

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
for ( s in sum_list ) {
source(snakemake@config[["differential_plots"]][["scripts"]][[paste(s)]])
}
# GC, Length and RPKM Bias ===================================
source(snakemake@config[["differential_plots"]][["scripts"]][["bias"]])


# Heatmap
if (difference != "splicing ratio") {
source(snakemake@config[["differential_plots"]][["scripts"]][["heat"]])
}

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

if (difference == "splicing ratio") {
plots <- list("summary","bias","ma_plot","volcano_plot")
} else {
plots <- list("summary","bias","ma_plot","volcano_plot","meta_plot","heatmap")
}

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
doc <- body_add(doc,value=plot_p,width = 6, height = 7.5, res= plot_dpi,style = "centered")
}

doc <- body_add(doc,title_p)
# Loop for each line in caption
for (c in caption_p) {
doc <- body_add(doc,fpar(values=c,fp_p = fp_par(padding.top=(caption_size/2))))
}

if ( (plot_n-1) == length(plots)) {
  doc <- body_add(doc,run_pagebreak())
}

}
}
}

sum_bar_data

head(sum_bar_data,5)
head(sum_pie_data,5)
head(sum_violin_data,5)

source(snakemake@config[["differential_plots"]][["scripts"]][["overview"]])
#doc <- cursor_begin(doc)
head(expr_i,10)
head(heat_data,10)

print(doc,target=snakemake@output[["docx"]])

doc <- read_docx(snakemake@input[["docx"]])
analysis_heading <- paste( "Differential", toTitleCase(feature), toTitleCase(difference))
doc <- body_add(doc,fpar(ftext(analysis_heading, prop=heading_2)),style = "heading 2")

doc <- body_add(doc,fpar(ftext("Overview", prop=heading_3)),style = "heading 3")
#doc <- body_add(doc,value=overview,width = 6, height = 9, res= plot_dpi,style = "centered")
doc <- body_add(doc,fpar(ftext("A", prop=bold)),style = "Normal")
doc <- body_add(doc,value=sum_pie,width = 6, height = 2.7, res= plot_dpi,style = "centered")

doc <- body_add(doc,fpar(ftext("B", prop=bold)),style = "Normal")
doc <- body_add(doc,value=sum_bar,width = 6, height = 2.7, res= plot_dpi,style = "centered")

doc <- body_add(doc,fpar(ftext("C", prop=bold)),style = "Normal")
doc <- body_add(doc,value=sum_violin,width = 6, height = 2.7, res= plot_dpi,style = "centered")

doc <- body_add(doc,run_pagebreak())
doc <- body_add(doc,block_pour_docx(snakemake@output[["docx"]]))

print(doc, target = snakemake@output[["docx"]])

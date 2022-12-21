log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(dplyr)
library(tools)
library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)
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

#sig_bg <- read.csv(snakemake@input[["sig_bg"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
sig <- as.numeric(snakemake@params[["sig"]])
bg <- as.numeric(snakemake@params[["bg"]])
section <- snakemake@params[["section"]]
base <- gsub("_"," ",gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@params[["base"]])))
base_feat <- gsub("_"," ",gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@params[["base_feat"]])))

mx_samples <- c(snakemake@params[["samples"]])

use_base_length <- (section!="body" & as.logical(snakemake@params[["main_int"]])) 
use_base_length

# Import snakemake parameters and inputs
# Identify analysis
diff <- as.character(snakemake@wildcards[["difference"]])
difference <- as.character(diff)

analysis <- paste("differential", difference, "analysis")

biases <- read.csv(snakemake@config[["bias"]],header=T, sep='\t', check.names=FALSE)
biases
paste(biases$bias,"_count_bias",sep='')

if (difference != "splicing_ratio") {

  is_antisense <- snakemake@params[["is_antisense"]]==-1
  difference <- gsub("[_.]"," ",diff)
  difference_unit <- paste(difference, "(RPKM)")
  sig_bg <- read.csv(snakemake@input[["sig_bg"]],header=T,row.names = 1, sep='\t', check.names=FALSE)

  measurement <- "Count"

  mx_df <- read.csv(snakemake@input[["mx_data"]],header=T, sep='\t', check.names=FALSE)
  bef_bin <- snakemake@params[["bef_bin"]]
  main_bin <- snakemake@params[["main_bin"]]
  section <- snakemake@params[["section"]]
  plot_median <- as.logical(snakemake@config[["metagene"]][["plot_median"]])
  meta_y <- ifelse(plot_median,"Median","Mean")
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

plotbef_brk_len <- signif(plotbef_len*1.5,1)/2
plotaft_brk_len <- signif(plotaft_len*1.5,1)/2
plotbef_brk_pos <- 0-bef_bin-(signif(plotbef_bin*1.5,1)/2)
plotaft_brk_pos <- main_bin-bef_bin + (signif(plotaft_bin*1.5,1)/2)
plotbef_brk <- paste(ifelse(is_antisense,"+","-"), as.character(plotbef_brk_len), sep="")
plotaft_brk <- paste(ifelse(is_antisense,"-","+"), as.character(plotaft_brk_len), sep="")

if (section=="body") {
meta_xbrks <- c(
 if(plotbef_len>0) (plotbef_brk_pos) else NULL, 
 0,
 main_bin, 
 if(plotaft_len>0) (plotaft_brk_pos) else NULL
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

} else {
  difference <- gsub("_"," ",difference)
  difference_unit <- gsub("_"," ",difference)
  measurement <- diff
}

# Import wildcards as text
prefix <- gsub("([[:lower:]])([[:upper:]])",perl=TRUE,"\\1 \\2",as.character(snakemake@wildcards[["prefix"]]))
tag <- toTitleCase(as.character(snakemake@wildcards[["tag"]]))
valid <- toTitleCase(as.character(snakemake@wildcards[["valid"]]))
feature <- gsub("_"," ",gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["feature"]])))
title_feature <- gsub("(?<!\\w)(.)","\\U\\1", feature, perl = TRUE)
title_base <- paste("Base",gsub("(?<!\\w)(.)","\\U\\1", base, perl = TRUE))
protocol <- as.character(snakemake@params[["protocol"]])

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

spikein <- gsub("_"," ",as.character(snakemake@wildcards[["spikein"]]))
spikein <- if (difference == "splicing ratio") ("internal") else (spikein)

contrast <- list(c(control,treat))


splice <- gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@wildcards[["splice"]]))
splice <- if (difference=="splicing ratio") ("All") else (splice)
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

genesets

gene_sets

goi <- unlist(lapply(unlist(strsplit(paste(gene_sets$genes[(gene_sets$set_name %in% genesets) & !as.logical(gene_sets$sep)],collapse=","),",")),trimws))
sep_gene_sets <- gene_sets$set_name[(gene_sets$set_name %in% genesets) & as.logical(gene_sets$sep)]

sep_gene_sets

analysis_title <- paste(
  common, ifelse(common=="",""," "), treat, " vs ", control, " ", protocol, " ", prefix, " ", ifelse(difference != "splicing ratio", paste(spikein, " ", normaliser, " normalised ",sep=""),""), difference,
#  common, ifelse(common=="",""," "), treat, " vs ", control, " ", protocol, "\n",
#  prefix, " ", ifelse(difference != "splicing ratio", paste(spikein, " ", normaliser, " normalised "),""), difference, "\n",
#  valid, " ", tag, " ", feature,
  sep="")

analysis_title

#goi <- c(snakemake@params[["GOI"]])

#  Import sample config, size factors and count table and conduct DESeq2 Differential Expression Analysis
rep_pair <- as.logical(snakemake@params[["paired"]])

# Import sample table and define ColData dataframe for deseq2
sample_table <- read.table(snakemake@config[["samples"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table$sample_name <- paste(sample_table$condition,"_",sample_table$protocol,"_Replicate_",sample_table$replicate,sep="")
sample_table <- sample_table[sample_table$condition %in% c(control_cond, treatment),]
sample_table <- sample_table[sample_table$protocol == as.character(snakemake@params[["protocol"]]),]
rownames(sample_table) <- sample_table$sample_name

mean_level <- read.csv(snakemake@input[["levels"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
names(mean_level) <- c(control,treat)
cts <- read.csv(snakemake@input[["counts"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
cts <- cts[names(cts) %in% sample_table$sample_name]

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
measurements <- unlist(lapply(rownames(analysis_plots),trimws))
analysis_plots <- data.frame(lapply(analysis_plots,trimws),check.names=F)
rownames(analysis_plots) <- measurements
measurement
sum_list <- as.list(colnames(analysis_plots)[analysis_plots[rownames(analysis_plots)==measurement,]==TRUE])
sum_list
sum_list <- c(sum_list,"ma","volcano")


# Initialise Word Document
doc <- read_docx(snakemake@input[["docx"]])
min_rpkm_pc <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_rpkm_percentile"]])

# Plot figures for features in each mono/multiexonic-biotype groups
i_group <- c(c(""),unique(expr$group[expr$biotype %in% biotypes]),sep_gene_sets)
for (i in i_group) {

#i <- "multiexonic protein coding"

#i <- "Histones"

if (i =="") {
  expr_i <- expr
} else if (i %in% gene_sets$set_name) {
  i_set_genes <- unlist(lapply(unlist(strsplit(paste(gene_sets$genes[gene_sets$set_name==i],collapse=","),",")),trimws))
  expr_i <- expr[expr$root_name %in% i_set_genes,]
} else if (i %in% expr$group) {
  expr_i <- expr[expr$group==i,]
} 

expr_full <- expr_i

i <- gsub("_"," ",i)

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

cts_i <- cts[rownames(cts) %in% expr_i$featureID,unlist(lapply(cts,is.numeric))] 

if (as.logical(snakemake@config[["differential_analysis"]][["use_p_adj_min_mean"]])) {
  min_mean <- max(expr_i$baseMean[is.na(expr_i$padj)])
} else {
  min_mean <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_mean_reads"]])
}

config_min_mean <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_mean_reads"]])

#meta_gene_n <- sum(rownames(sig_bg)[ (sig_bg$sig2bg >= sig & sig_bg$bg2sig >= bg) ] %in% expr_i$featureID[expr_i$baseMean >= snakemake@config[["metagene"]][["min_reads"]]])

#heat_min_rpkm <- quantile(expr_i$RPKM[expr_i$baseMean > 0],snakemake@config[["heatmap"]][["min_rpkm_pc"]]/100,na.rm=T)
heat_min_reads <- snakemake@config[["heatmap"]][["min_reads"]] 
expr_heat <- expr_i[expr_i$baseMean >= heat_min_reads,]
expr_heat$padj <- ifelse(is.na(expr_heat$padj),expr_heat$pvalue,expr_heat$padj)
expr_heat$padj <- ifelse(is.na(expr_heat$padj),1,expr_heat$padj)
expr_heat <- expr_heat %>% mutate(p_rank=1) %>% ungroup
expr_heat <- expr_heat %>% mutate(lfc_rank=1) %>% ungroup
expr_heat$log10P <- -log10(expr_heat$padj)
expr_heat$colour <- ifelse(expr_heat$padj < sig_p, ifelse(expr_heat$log2FoldChange < 0, down_col, up_col), insig_col)


total_i <- nrow(expr_heat)

min_rpkm <- quantile(expr_i$RPKM[expr_i$baseMean > 0],min_rpkm_pc/100,na.rm=T)
#meta_gene_n <- sum(rownames(sig_bg)[ (sig_bg$sig2bg >= sig & sig_bg$bg2sig >= bg) ] %in% expr_i$featureID[expr_i$baseMean >= config_min_mean & expr_i$baseMean >= min_rpkm])

#insuf_i <- sum( (expr_heat$baseMean < min_mean) | (!is.na(expr_heat$RPKM) & (expr_heat$RPKM < min_rpkm)) | (is.na(expr_heat$pvalue)),na.rm=T )

expr_i <- expr_i[(expr_i$baseMean >= min_mean & expr_i$RPKM >= min_rpkm),]


expr_i$padj <- ifelse(is.na(expr_i$padj),expr_i$pvalue,expr_i$padj)
expr_i <- expr_i[!is.na(expr_i$padj),]

feature_i <- ifelse(i=="", feature, paste(i,feature))
base_i <- ifelse(i=="", paste("base", base), paste("base", i,base)) 

title_feature_i <- gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE)
title_base_i <- gsub("(?<!\\w)(.)","\\U\\1", base_i, perl = TRUE)

# Write heading for this analysis group
group_heading <- toTitleCase(ifelse(i=="","All",i))
doc <- body_add(doc,fpar(ftext(group_heading, prop=heading_3)),style = "heading 3")

group_label <- gsub("exonic ","exonic\n",group_heading)

head(expr_i,5)
if (nrow(expr_i)==0 | nrow(expr_heat)==0 ) {

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

lfc_max <- max(c(quantile(abs(expr_i$log2FoldChange),lfc_pc_lim/100,na.rm=T),0))
expr_i$colour <- ifelse(expr_i$padj < sig_p, ifelse(expr_i$log2FoldChange < 0, down_col, up_col), insig_col)

title_i <- gsub("_"," ",paste(experiment, feature_i, "Differential", toTitleCase(difference),sep=" "))
head(expr_heat,10)
head(expr_i,10)

# Summary Plots =============================================
for ( s in sum_list ) {
source(snakemake@config[["differential_plots"]][["scripts"]][[paste(s)]])
}
# GC, Length and RPKM Bias ===================================
source(snakemake@config[["differential_plots"]][["scripts"]][["bias"]])
head(expr_bias,10)
source(snakemake@config[["differential_plots"]][["scripts"]][["count_bias"]])


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
summary_h <- 8

plots <- c("summary","bias","count_bias","ma_plot","volcano_plot")

if (difference != "splicing ratio") {
plots <- c(plots,"meta_plot","heatmap")
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
doc <- body_add(doc,fpar(ftext(toTitleCase(gsub("_"," ",p)), prop=heading_3)),style = "heading 4")
doc <- body_add(doc,value=plot_p,width = 6, height = get(paste(p,"_h",sep="")), res= plot_dpi,style = "centered")
if (p=="summary") {
doc <- body_add(doc,run_pagebreak())
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


#=====Test====

#=====testEnd======


count_bias
head(count_bias_data,5)
count_bias_title
count_bias_h
count_bias_caption
sum_bar_data

head(sum_bar_data,5)
sum_pie_data
head(sum_violin_data,5)
write.table(sum_bar_data,file=snakemake@params[["bar_data"]],sep='\t',row.names=F,quote=F)
write.table(sum_pie_data,file=snakemake@params[["pie_data"]], sep='\t',row.names=F,quote=F)
write.table(sum_violin_data,file=snakemake@params[["violin_data"]],sep='\t',row.names=F,quote=F)

source(snakemake@config[["differential_plots"]][["scripts"]][["overview"]])
#doc <- cursor_begin(doc)
head(expr_i,10)
head(heat_data,10)

i_group_levels

print(doc,target=snakemake@output[["docx"]])

doc <- read_docx(snakemake@input[["docx"]])
analysis_heading <- paste( "Differential", toTitleCase(feature), toTitleCase(difference))
doc <- body_add(doc,fpar(ftext(analysis_heading, prop=heading_2)),style = "heading 2")

doc <- body_add(doc,fpar(ftext("Overview", prop=heading_3)),style = "heading 3")
#doc <- body_add(doc,value=overview,width = 6, height = 9, res= plot_dpi,style = "centered")
doc <- body_add(doc,fpar(ftext("A", prop=bold)),style = "Normal")
doc <- body_add(doc,value=sum_pie,width = 6, height = 2.2, res= plot_dpi,style = "centered")

doc <- body_add(doc,fpar(ftext("B", prop=bold)),style = "Normal")
doc <- body_add(doc,value=sum_bar,width = 6, height = 2.2, res= plot_dpi,style = "centered")

doc <- body_add(doc,fpar(ftext("C", prop=bold)),style = "Normal")
doc <- body_add(doc,value=sum_violin,width = 6, height = 2.2, res= plot_dpi,style = "centered")

#doc <- body_add(doc,run_pagebreak())

doc <- body_add(doc,block_pour_docx(snakemake@output[["docx"]]))

doc <- headers_replace_all_text(doc,"",analysis_title)

print(doc, target = snakemake@output[["docx"]])
save.image(file = snakemake@output[["rdata"]])

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(dplyr)
library(tools)
library(tidyverse)
library(ggrepel)
library(ggtext)
library(stringr)
library(rvg)
library(scales)


parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Initialise experiment and feature Settings
expr <- read.csv(snakemake@input[["lfc"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
expr$featureID <- rownames(expr)
mean_level <- read.csv(snakemake@input[["levels"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
cts <- read.csv(snakemake@input[["counts"]],header=T,row.names = 1, sep='\t', check.names=FALSE)

genes <- read.csv(snakemake@input[["genetab"]],sep='\t',header=F,col.names=c("gene_id","gene_name","biotype","exon_count","chr","start","end","strand"),check.names=F)
bed <- read.csv(snakemake@input[["bed"]],header=F, sep='\t', check.names=FALSE)[,c(4,8,11)]
names(bed) <- c("rootID","featureID","baseID")
expr$baseID <- bed$baseID[match(expr$featureID,bed$featureID)]
expr$rootID <- bed$rootID[match(expr$featureID,bed$featureID)]
expr$root_name <- genes$gene_name[match(expr$rootID,genes$gene_id)]
head(bed,10)

sense_ls <- snakemake@input[["sense_mx"]]
antisense_ls <- snakemake@input[["antisense_mx"]]
sig_bg <- read.csv(snakemake@input[["sig_bg"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
sig <- as.numeric(snakemake@params[["sig"]])
bg <- as.numeric(snakemake@params[["bg"]])

if (length(antisense_ls)>0) {
sense_dirs <- c("sense","antisense")
} else {
sense_dirs <- c("sense")
}


plot_median <- as.logical(snakemake@config[["coverage_plots"]][["metagene"]][["plot_median"]])

head(sig_bg,10)
meta_features <- rownames(sig_bg)[ (sig_bg$sig2bg >= sig & sig_bg$bg2sig >= bg) ] 
head(meta_features,10)
mx_samples <- c(snakemake@params[["samples"]])

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

# Import wildcards as text
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

#goi <- c(snakemake@params[["GOI"]])

# Import sample table and define ColData dataframe for deseq2
sample_table <- read.table(snakemake@config[["samples"]], sep='\t',header=TRUE, check.names=FALSE)
sample_table$sample_name <- paste(sample_table$condition,"_",sample_table$protocol,"_Replicate_",sample_table$replicate,sep="")
sample_table <- sample_table[sample_table$condition %in% c(control_cond, treatment),]
sample_table <- sample_table[sample_table$protocol == as.character(snakemake@params[["protocol"]]),]
rownames(sample_table) <- sample_table$sample_name

sample_table$cond_abbr <- c(control,treat)[match(sample_table$condition, c(control_cond,treatment))]
sample_table$sample_abbr <- paste(sample_table$cond_abbr," Rep",sample_table$replicate,sep="") 


# Initialise Plotting
for (i in unique(expr$biotype[!is.na(expr$exon)])) {
  print(paste(i, unique(expr$exon[!is.na(expr$exon) & expr$biotype==i]),length(unique(expr$exon[!is.na(expr$exon) & expr$biotype==i]))))
  if (length(unique(expr$exon[!is.na(expr$exon) & expr$biotype==i]))==1) {
    expr$group[expr$biotype==i] <- gsub("_", " ", i)
  }
}

min_rpkm_pc <- as.numeric(snakemake@config[["differential_analysis"]][["minimum_rpkm_percentile"]])
# Import size factors
size_table <- read.csv(snakemake@input[["size_table"]],header=T,sep="\t",check.names=F)
size_table <- size_table[match(mx_samples,size_table$sample_name),]

# Import metaplot params and matrices
plotbef_bin <- snakemake@params[["plotbef_bin"]]
plotaft_bin <- snakemake@params[["plotaft_bin"]]
base <- as.character(snakemake@params[["base"]])
bef_bin <- snakemake@params[["bef_bin"]]
main_bin <- snakemake@params[["main_bin"]]
section <- snakemake@params[["section"]]
len_bef_n <- as.numeric(snakemake@params[["len_bef"]])
len_aft_n <- as.numeric(snakemake@params[["len_aft"]])
len_bef <- paste("-",as.character(snakemake@params[["len_bef"]]),sep="")
len_aft <- paste("+",as.character(snakemake@params[["len_aft"]]),sep="")

for ( i in sense_dirs ) {
mx_ls <- get(paste(i,"_ls",sep=""))
mx <- lapply(mx_ls,function(x) {read.table(x,sep='\t',header=T,row.names=1,check.names=FALSE)})
mx <- lapply(mx,function(x) {x[rownames(x) %in% meta_features,]})
mx <- lapply(1:length(mx),function(x) {as.data.frame(mx[[x]])*size_table$scale_factor[size_table$sample_name==mx_samples[x]]
})
names(mx) <- sample_table$sample_abbr[match(mx_samples,sample_table$sample_name)]
assign(paste(i,"_mx",sep=""),mx)
}

head(sense_mx[[1]][,1:10],10)

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

i <- gsub("_"," ",i)

head(expr_i, 5)

min_mean <- as.numeric(snakemake@config[["coverage_plots"]][["metagene"]][["min_reads"]])

min_rpkm <- if (min_rpkm_pc==0) (0) else (as.numeric(quantile(expr_i$RPKM[expr_i$baseMean > 0],min_rpkm_pc/100)))
min_mean
min_rpkm


expr_i <- expr_i[( (expr_i$baseMean >= min_mean) & (expr_i$RPKM >= min_rpkm) ),]

#expr_i$padj <- ifelse(is.na(expr_i$padj),expr_i$pvalue,expr_i$padj)
#expr_i <- expr_i[!is.na(expr_i$padj),]

head(expr_i,10)

meta_trim <- as.numeric(snakemake@config[["coverage_plots"]][["metagene"]][["anomaly_trim"]])

# Process matrices to trimmed mean and/or/maybe spline
for ( s in sense_dirs ) {

#s <- "sense"

mx <- get(paste(s,"_mx",sep=""))
head(mx[[1]][,1:10],10)

mx_names <- rownames(mx[[1]])[1:20]
mx_names
expr_names <- rownames(expr_i)[1:20]
expr_names
match(mx_names,expr_names)
ls <- (rownames(mx[[1]]) %in% rownames(expr_i))
ls[1:20]
mx <- lapply(mx,function(x) {x[match(expr_i$featureID,rownames(x),]} )
head(mx[[1]][,1:10],10)

# Plot Coverage Heatmaps, against rankings defined in heat_config


# Coverage Heatmaps
mean_rpkm <- mean(expr_i$RPKM, na.rm=T)
rpkm_norm_factor <- if (mean_rpkm > 0) mean_rpkm/expr_i$RPKM
norm_mx <- lapply(1:length(mx),function(x) {as.data.frame(mx[[x]])*rpkm_norm_factor})
names(norm_mx) <- sample_table$sample_abbr[match(mx_samples,sample_table$sample_name)]

#mx_sum_median <- data.frame(lapply(mx, function(x) {apply(x,2,function(y) {median(y,na.rm=F)})}),check.names=F)
#mean_mx <- data.frame(lapply(mx, function(x) {apply(x,2,function(y) {mean(y,na.rm=F,trim=meta_trim)})}),check.names=F)
sum_mx <- data.frame(lapply(mx, function(x) {apply(x,2,function(y) {sum(y,na.rm=T)})}),check.names=F)

head(sum_mx,10)

pos <- as.numeric(rownames(sum_mx))
rev_pos <- rev(as.numeric(rownames(sum_mx)))

if (s=="sense") {
pos_s <- pos
} else {
pos_s <- rev_pos
}

pos_s

mx_data <- data.frame(
  unlist(
    lapply(colnames(sum_mx), function (x) {
      ifelse(s=="sense", sum_mx[paste(x)], 0-sum_mx[paste(x)])
    })),
  unlist(
    lapply(colnames(sum_mx), function (x) {
      replicate(nrow(sum_mx),paste(x))
    })),
  c(rep(pos_s,ncol(sum_mx))),
  unlist(
    lapply(colnames(sum_mx), function (x) {
      replicate(nrow(sum_mx),paste(x,s,sep="_"))
    }))
)
names(mx_data) <- c("value","Sample","Position","variable")

mx_data$sense <- s
assign(paste(s,"_mx_data",sep=""),mx_data)
}

head(sense_mx_data,10)

if (length(antisense_ls)>0) {
sense_mx_data <- rbind(sense_mx_data,antisense_mx_data)
}

sense_mx_data$Condition <- sample_table$cond_abbr[match(sense_mx_data$Sample,sample_table$sample_abbr)]
sense_mx_data$cond_group <- paste(sense_mx_data$Condition, sense_mx_data$sense)
sense_mx_data$colour <- sample_table$colour[match(sense_mx_data$Sample,sample_table$sample_abbr)]
sample_colours <- as.character(unique(sense_mx_data$colour))
sample_names <- as.character(unique(sense_mx_data$Sample))
sample_colours
sample_names
names(sample_colours) <- gsub("_"," ",sample_names)
sample_colours

sense_mx_data$rep <- paste("Replicate", sample_table$replicate[match(sense_mx_data$Sample,sample_table$sample_abbr)])
sense_mx_data$Sample <- gsub("_"," ",sense_mx_data$Sample)
sense_mx_data$i_group <- i 
head(sense_mx_data,10)



if (exists("mx_dataframe")) {
mx_dataframe <- rbind(mx_dataframe,sense_mx_data)
} else {
mx_dataframe <- sense_mx_data

}

}


write.table(mx_dataframe,file=snakemake@output[["mx_data"]],sep='\t',row.names=F,quote=F)
save.image(file = snakemake@output[["rdata"]])



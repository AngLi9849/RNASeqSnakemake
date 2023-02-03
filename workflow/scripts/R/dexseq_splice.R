log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(error=function()traceback(2))

library(DEXSeq)
library(ashr)
library(dplyr)
library(tools)

save_rdata <- as.logical(snakemake@config[["differential_analysis"]][["save_rdata"]])
parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Set Control and Treatment Conditions
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

# Import feature lengths and nucleotide content
express <- read.table(snakemake@input[["express"]], sep='\t',header=TRUE, check.names=FALSE)
nuc <- read.csv(snakemake@input[["nuc"]],header=T,row.names = 1, sep='\t')


# Import feature biotype and exon count information
genes <- read.csv(snakemake@input[["genetab"]],sep='\t',header=F,col.names=c("geneID","gene_name","biotype","exon_count","chr","start","end","strand"),check.names=F)
bed <- unique(read.csv(snakemake@input[["bed"]],sep='\t',header=F,check.names=F)[,c(4,8,9,11)])
names(bed) <- c("rootID","featureID","feature_name","baseID")

genes <- genes[genes$geneID %in% bed$rootID,]
genes[,c("rootID","featureID","feature_name","baseID")] <- bed[match(genes$geneID,bed$rootID),c("rootID","featureID","feature_name","baseID")]

rep_pair <- as.logical(snakemake@params[["paired"]])

#Import splice and unsplice counts
splice_cts <- read.table(snakemake@input[["spliced"]], header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)
unsplice_cts <- read.table(snakemake@input[["unspliced"]], header=TRUE, row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)

splice_cts <- splice_cts[ , order(names(splice_cts))]
unsplice_cts <- unsplice_cts[ , order(names(unsplice_cts))]
samples <- names(splice_cts)

# Import sample table
sample_table <- read.table(snakemake@config[["samples"]], sep='\t',header=TRUE, check.names=FALSE,comment.char="")
sample_table$sample_name <- paste(sample_table$condition,"_",sample_table$protocol,"_Replicate_",sample_table$replicate,sep="")
sample_table <- sample_table[sample_table$sample_name %in% samples,]
rownames(sample_table) <- sample_table$sample_name
sample_table <- sample_table[order(row.names(sample_table)), , drop=F]

splice_cts <- splice_cts[,names(splice_cts) %in% sample_table$sample_name]
unsplice_cts <- unsplice_cts[,names(unsplice_cts) %in% sample_table$sample_name]

sample_table <- sample_table[match(names(splice_cts),sample_table$sample_name),]

#Count table of unnormalised spliced and spliced counts
cts <- rbind(splice_cts,unsplice_cts)
cts$state <- c(replicate(nrow(splice_cts),"spliced"),replicate(nrow(unsplice_cts),"unspliced"))
cts$id <- c(rownames(splice_cts),rownames(splice_cts))

splice_cts_names <- row.names(splice_cts)
splice_cts <- sapply(splice_cts,as.numeric)
row.names(splice_cts) <- splice_cts_names

unsplice_cts_names <- row.names(unsplice_cts)
unsplice_cts <- sapply(unsplice_cts,as.numeric)
row.names(unsplice_cts) <- unsplice_cts_names

#Internally Normalised per-gene splice and unspliced counts 
cts_sum <- splice_cts + unsplice_cts
cts_mean <- apply(cts_sum,1,mean)
cts_mean <- ifelse(is.na(splice_cts),cts_mean,cts_mean)
unsplice_cts <- round(ifelse((cts_mean<=1 | cts_sum <= 1),unsplice_cts,unsplice_cts*cts_mean/cts_sum))
splice_cts <- round(ifelse((cts_mean<=1 | cts_sum <= 1),splice_cts,splice_cts*cts_mean/cts_sum))


#splice_mean <- apply(splice_cts,1,mean)
#splice_mean <- ifelse(is.na(splice_cts),splice_mean,splice_mean)
#unsplice_cts <- round(ifelse((splice_mean<=1 | splice_cts <= 1),unsplice_cts,unsplice_cts*splice_mean/splice_cts))
#splice_cts <- round(ifelse((splice_mean<=1 | splice_cts <= 1),splice_cts,splice_mean))

feature_id <- row.names(splice_cts)
group_id <- rep("input",length(feature_id))

#feature_id <- c(row.names(splice_cts),paste(row.names(unsplice_cts),"_unsplice",sep=""))
#group_id <- rep(row.names(splice_cts),2)

coldata <- sample_table[,c("condition","replicate")]

splice_control <- splice_cts[,sample_table$condition==control_cond]
unsplice_control <- unsplice_cts[,sample_table$condition==control_cond]
coldata_control <- coldata[sample_table$condition==control_cond,]

# Loop for each expriment condition against control
#for (a in experiments) {

splice_exp <- splice_cts[,sample_table$condition==treatment]
unsplice_exp <- unsplice_cts[,sample_table$condition==treatment]

splice_cts_exp <- cbind(splice_exp,splice_control)
unsplice_cts_exp <- cbind(unsplice_exp,unsplice_control)

#row.names(unsplice_cts_exp) <- paste(row.names(unsplice_cts_exp),"_unsplice",sep="")
#cts_exp <- rbind(splice_cts_exp,unsplice_cts_exp)
#head(cts_exp,10)

coldata_exp <- coldata[sample_table$condition==treatment,]

coldata <- rbind(coldata_exp,coldata_control)
coldata_samples <- rownames(coldata)

coldata <- data.frame(lapply(coldata,function(x) { gsub("[\\+|-|_]",".",x) } ))
rownames(coldata) <- coldata_samples
coldata$condition <- factor(coldata$condition, levels=c(gsub("[\\+|-|_]",".",control_cond),gsub("[\\+|-|_]",".",treatment)))

control_cond
treatment
sample_table
coldata
head(splice_cts_exp,10)
nrow(splice_cts_exp)
head(unsplice_cts_exp,10)
nrow(unsplice_cts_exp)


if (rep_pair){
  full_model <- ~sample + exon + condition:exon + replicate:exon
  reduced_model <- ~sample + exon + replicate:exon
 } else {
  full_model <- ~sample + exon + condition:exon
  reduced_model <- ~sample + exon
}

dds <- DEXSeqDataSet(countData=splice_cts_exp,sampleData=coldata, featureID=feature_id, groupID=group_id, alternativeCountData = unsplice_cts_exp, design=full_model)
#dds <- DEXSeqDataSet(countData=cts_exp,sampleData=coldata, featureID=feature_id, groupID=group_id, design=full_model)

rds <- estimateSizeFactors(dds)
rds <- estimateDispersions(rds,formula=full_model)
rds <- testForDEU(rds,reducedModel = reduced_model,fullModel = full_model)
rds <- estimateExonFoldChanges(rds, fitExpToVar="condition")
res <- DEXSeqResults(rds)

expr <- data.frame(res@listData[1:10],check.names=FALSE)
#expr <- expr[unlist(lapply(expr$featureID,function(x) {!grepl(x,"_unsplice")})),]
rownames(expr) <- splice_cts_names 
colnames(expr)[10] <- "Rawlog2FoldChange"
colnames(expr)[3] <- "baseMean"
#expr <- expr[!is.na(expr$padj),]
expr[,c("featureID","feature_name","rootID","root_name","baseID","exon_count","biotype")] <- genes[match(rownames(expr),genes$featureID),c("featureID","feature_name","geneID","gene_name","baseID","exon_count","biotype")]
expr$exon <- ifelse(expr$exon_count>1,"multiexonic","monoexonic")
expr$change <- ifelse(expr$Rawlog2FoldChange>=0,"Increased","Decreased")
expr$group <- gsub("_"," ",paste(expr$exon,expr$biotype))
expr$group2<-paste(expr$change,expr$group)
expr$log2FoldChange <- expr$Rawlog2FoldChange
expr$featureID <- rownames(expr)

#Use base feature length if lengths of tested feature is fixed
section <- snakemake@params[["section"]]
base_bed <- read.csv(snakemake@input[["base_bed"]],header=F, sep='\t', check.names=FALSE)[,c(5,8)]
names(base_bed) <- c("Length","baseID")
#base_bed <- data.frame("baseID" = unique(base_bed$baseID), "Length"=lapply(unique(base_bed$baseID) function(x) {sum(base_bed$Length[base_bed$baseID==x])}), check.names=F)
use_base_length <- (section!="body" & as.logical(snakemake@params[["main_int"]]))
head(base_bed,10)
use_base_length
if (use_base_length) {
expr$Length <- base_bed$Length[match(expr$baseID,base_bed$baseID)]
} else {
expr$Length <- express$Length[match(rownames(expr),express$featureID)]
}

expr[,c("GC","AT")] <- nuc[match(rownames(expr),rownames(nuc)),c("GC","AT")]
expr$RPKM <- express$RPKM[match(rownames(expr),express$featureID)]
expr$RPK <- express$RPK[match(rownames(expr),express$featureID)]



#expr$log10P <- -log10(expr$padj)
#expr$z <- qnorm(expr$pvalue)
#expr$lfcSE <- abs(expr$Rawlog2FoldChange/expr$z)
#ash <- ash(betahat=expr$Rawlog2FoldChange,sebetahat=expr$lfcSE, mixcompdist = "normal",method= "shrink")
#expr$log2FoldChange <- ash$result$PosteriorMean
#expr <- expr %>% arrange(group2,padj) %>% group_by(group2) %>% mutate(p_rank=1:n()) %>% ungroup
#expr <- expr %>% arrange(group2,abs(log2FoldChange)) %>% group_by(group2) %>% mutate(lfc_rank=n():1) %>% ungroup

#expr <- expr %>% arrange(change,padj) %>% group_by(change) %>% mutate(total_p_rank=1:n()) %>% ungroup
#expr <- expr %>% arrange(change,abs(log2FoldChange)) %>% group_by(change) %>% mutate(total_lfc_rank=n():1) %>% ungroup

# Total Splicing Ratio and mean feature splicing ratios
splice_ratio <- data.frame(splice_cts/(splice_cts + 2*unsplice_cts),check.names=F)
head(splice_ratio,10)

splice_ratio_mean <- data.frame(lapply(c(control_cond,treatment),function(x) {
    if ( sum(colnames(splice_ratio) %in% sample_table$sample_name[sample_table$condition==x])==1 ) {
      splice_ratio[,colnames(splice_ratio) %in% sample_table$sample_name[sample_table$condition==x]] 
    } else {
      apply(splice_ratio[,colnames(splice_ratio) %in% sample_table$sample_name[sample_table$condition==x]],1,FUN=mean)
    }
  }))
names(splice_ratio_mean) <- c(control_cond,treatment)
splice_ratio_mean <- splice_ratio_mean[rownames(splice_ratio_mean) %in% expr$featureID,]
head(splice_ratio_mean,10)

mean_level <- splice_ratio_mean[complete.cases(splice_ratio_mean),]

write.table(data.frame("id"=rownames(expr),expr, check.names=FALSE),file=snakemake@output[["lfc"]],sep='\t',row.names=F,quote=F)

write.table(data.frame("id"=rownames(mean_level),mean_level, check.names=FALSE),file=snakemake@output[["levels"]],sep='\t',row.names=F,quote=F)

write.table(data.frame("id"=rownames(cts),cts, check.names=FALSE),file=snakemake@output[["counts"]],sep='\t',row.names=F,quote=F)


#}
if (save_rdata) {
  save.image(file = snakemake@params[["rdata"]])
}

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log,type="output")
sink(log,type="message")

library("Rcpp")
library("tidyverse")
library("reshape2")
library("tools")

# Import scaled matrices and sample names as lists and samples config table
sample_table <- read.table(snakemake@input[["sample_table"]],header=T,sep='\t')
matrices_list <- snakemake@input[["matrices"]]
matrices_df <- lapply(matrices_list,function(x) {read.table(x,sep='\t',header=F,skip=1)})
sample_names <- lapply(matrices_list, function(x) {readLines(x,n=1)})

# Set variables for output
outpdf <- snakemake@output[["meta_profile_pdf"]]
outpng <- snakemake@output[["meta_profile_png"]]

# Import matrix parameters as variables
feature <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["feature"]])))
score <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["score"]])))
bin_size <- as.numeric(snakemake@params[["bin_size"]])
before <- as.numeric(snakemake@params[["before"]])
start <- as.numeric(snakemake@params[["start"]])
body_length <- as.numeric(snakemake@params[["body_length"]])
end <- as.numeric(snakemake@params[["end"]])
after <- as.numeric(snakemake@params[["after"]])

# Normalise to trimmed size factor, calculate trimmed_means of each coordinate on each matrix and combine into a melted long table of metadata
trimmed_colmean <- lapply(matrices_df, function(x) {apply(x,2,function(y) {mean(y,na.rm=F,trim=snakemake@params[["trim"]])})})
trimmed_size <- lapply(trimmed_colmean,mean)
for(i in 1:length(trimmed_size)){trimmed_mean[[i]] <- trimmed_colmean[[i]]/trimmed_size[[i]]}
trimmed_mean <- data.frame(trimmed_mean)
colnames(trimmed_mean) <- sample_names
trimmed_mean$Position <- c(1:nrow(trimmed_mean)-(before+1))
meta_mean <- melt(trimmed_mean,id.vars="Position")
colnames(meta_mean) <- c("Position","Sample","value")

# Set meta-profile colours based on sample config table
meta_mean$colour <- sample_table$colour[match(meta_mean$Sample,sample_table$sample_name)]
sample_colours <- as.character(unique(meta_mean$colour))
sample_names <- as.character(unique(meta_mean$Sample))
names(sample_colours) <- sample_names

# Generate profile plot parameters
max_val <- ifelse(max(meta_mean$value) >= 0, signif(1.1*max(meta_mean$value),2), 0)
min_val <- ifelse(min(meta_mean$value) >= 0, 0 , signif(1.1*min(meta_mean$value),2))
xlim <- c((0-before),(start+body_length+end+after))
xbrks <- c(0,body_length+start+end,if(before>0)(signif(0-(before*1.5),1)/2) else NULL,if(after>0)(body_length+start+end+signif(after*1.5,1)/2) else NULL,if(start>0)start else NULL,if(end>0)(body_length+start) else NULL)
names(xbrks) <- c("TSS","TTS",if(before>0)(paste("-",signif(before*1.5,1)/2,sep="")) else NULL,if(after>0)(paste("+",signif(after*1.5,1)/2,sep="")) else NULL,if(start>0)(paste("+",start,sep="")) else NULL,if(end>0)(paste("-",end,sep="")) else NULL)
gene_number <- paste("n=",format(nrow(matrices_df[[1]]),big.mark=",",scientific=F),sep="")
title <- paste(feature," (",gene_number,") Normalised ",score,sep="")

# Plot meta-profile
meta_plot <- ggplot(meta_mean,mapping=aes(x=Position,y=value,colour=Sample)) + 
geom_line(size=0.5,stat="summary_bin",binwidth=bin_size) + 
scale_colour_manual(breaks=sample_names,values=sample_colours) + 
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) + 
geom_vline(xintercept=(body_length+start+end),color="black",alpha=0.2,linetype=5) + 
geom_vline(xintercept=start,color="black",alpha=0.2,linetype=5) +
geom_vline(xintercept=(body_length+start),color="black",alpha=0.2,linetype=5) +
scale_y_continuous(limits=c(min_val,max_val),breaks=c(min_val,min_val/2,0,max_val/2,max_val)) + 
geom_hline(yintercept=0,alpha=0.6) + 
scale_x_continuous(limits=xlim,breaks=xbrks) + 
ggtitle(title) + 
xlab("Relative Position (bps)") + 
ylab("Normalised Counts") + 
theme(panel.background=element_rect(fill="White",colour="white"), panel.border=element_rect(fill=NA,colour="black",size=0.7), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.7),axis.line.x.top=element_line(colour="black",size=0.7),axis.line.y.right=element_line(colour="black",size=0.7))
ggsave(filename=outpdf,plot=meta_plot,height=aspect*8,width=8,dpi=1200)
ggsave(filename=outpng,plot=meta_plot,height=aspect*8,width=8,dpi=1200)

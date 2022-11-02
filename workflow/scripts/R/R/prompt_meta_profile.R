log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log,type="output")
sink(log,type="message")

library("Rcpp")
library("tidyverse")
library("reshape2")

# Import scaled matrices and sample names as lists and samples config table
sample_table <- read.table(snakemake@input[["sample_table"]],header=T,sep='\t')

sense_matrices_list <- snakemake@input[["sense_matrices"]]
sense_matrices_df <- lapply(sense_matrices_list,function(x) {read.table(x,sep='\t',header=F,skip=1)})
sense_names <- lapply(sense_matrices_list, function(x) {readLines(x,n=1)})

antisense_matrices_list <- snakemake@input[["antisense_matrices"]]
antisense_matrices_df <- lapply(antisense_matrices_list,function(x) {read.table(x,sep='\t',header=F,skip=1)})
antisense_names <- lapply(antisense_matrices_list, function(x) {readLines(x,n=1)})

# Set variables for output
outpdf <- snakemake@output[["meta_profile_pdf"]]
outpng <- snakemake@output[["meta_profile_png"]]

# Import matrix parameters as variables
feature <- str_to_title(gsub("_"," ",as.character(snakemake@wildcards[["feature"]])))
bin_size <- as.numeric(snakemake@params[["bin_size"]])
prompt <- as.numeric(snakemake@params[["prompt"]])
tss <- as.numeric(snakemake@params[["tss"]])

# Normalise to sense strand total size, calculate trimmed means of each coordinate on each matrix and combine into a melted long table of metadata
sense_trimmed_colmean <- lapply(sense_matrices_norm, function(x) {apply(x,2,function(y) {mean(y,na.rm=F,trim=snakemake@params[["trim"]])})})
sense_trimmed_size <- lapply(sense_trimmed_colmean, function(x) {mean(x[1:(prompt-200)])})

sense_trimmed_mean <- NULL
for(i in 1:length(sense_trimmed_size)){sense_trimmed_mean[[i]] <- sense_trimmed_colmean[[i]]/sense_trimmed_size[[i]]}
sense_trimmed_mean <- data.frame(sense_trimmed_mean)
colnames(sense_trimmed_mean) <- sense_names
sense_trimmed_mean$Position <- c(1:nrow(sense_trimmed_mean)-(prompt+1))
sense_meta_mean <- melt(sense_trimmed_mean,id.vars="Position")
sense_meta_mean$Sample <- gsub("_sense","",sense_meta_mean$variable)

antisense_trimmed_colmean <- lapply(antisense_matrices_df, function(x) {apply(x,2,function(y) {mean(y,na.rm=F,trim=snakemake@params[["trim"]])})})
antisense_trimmed_mean <- NULL
for(i in 1:length(sense_trimmed_size)){antisense_trimmed_mean[[i]] <- antisense_trimmed_colmean[[i]]/sense_trimmed_size[[i]]}
antisense_trimmed_mean <- data.frame(antisense_trimmed_mean)
colnames(antisense_trimmed_mean) <- antisense_names
antisense_trimmed_mean$Position <- c(1:nrow(antisense_trimmed_mean)-(prompt+1))
antisense_meta_mean <- melt(antisense_trimmed_mean,id.vars="Position")
antisense_meta_mean$value <- 0-antisense_meta_mean$value
antisense_meta_mean$Sample <- gsub("_antisense","",antisense_meta_mean$variable)

meta_mean <- rbind(sense_meta_mean,antisense_meta_mean)
head(meta_mean,10)

# Set meta-profile colours based on sample config table
meta_mean$colour <- sample_table$colour[match(meta_mean$Sample,sample_table$sample_name)]
sample_colours <- as.character(unique(meta_mean$colour))
sample_names <- as.character(unique(meta_mean$Sample))
sample_colours
sample_names
names(sample_colours) <- sample_names
sample_colours

# Generate profile plot parameters
max_val <- ifelse(max(meta_mean$value) >= 0, signif(1.1*max(meta_mean$value),2), 0)
min_val <- ifelse(min(meta_mean$value) >= 0, 0 , signif(1.1*min(meta_mean$value),2))
xlim <- c((0-prompt),(prompt+tss))
xbrks <- c(0,if(prompt>0)(signif(0-(prompt*2),1)/2) else NULL,if(tss>0)(signif(tss*2,1)/2) else NULL)
names(xbrks) <- c("TSS",if(prompt>0)(paste("-",signif(prompt*2,1)/2,sep="")) else NULL,if(tss>0)(paste("+",signif(tss*2,1)/2,sep="")) else NULL)
gene_number <- paste("n=",format(nrow(sense_matrices_df[[1]]),big.mark=",",scientific=F),sep="")
title <- paste(feature," prompt-TSS (",gene_number,")",sep="")

# Plot meta-profile
meta_plot <- ggplot(meta_mean,mapping=aes(x=Position,y=value,group=variable,colour=Sample)) + 
geom_line(size=0.5,stat="summary_bin",binwidth=bin_size) + 
scale_colour_manual(breaks=sample_names,values=sample_colours) + 
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) + 
scale_y_continuous(limits=c(min_val,max_val),breaks=c(min_val,min_val/2,0,max_val/2,max_val)) + 
geom_hline(yintercept=0,alpha=0.6) + 
scale_x_continuous(limits=xlim,breaks=xbrks) + 
ggtitle(title) + 
xlab("Relative Position (bps)") + 
ylab("Normalised Counts") + 
theme(panel.background=element_rect(fill="White",colour="white"),panel.border=element_rect(fill=NA,colour="black",size=0.7), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.7),axis.line.x.top=element_line(colour="black",size=0.7),axis.line.y.right=element_line(colour="black",size=0.7))
ggsave(filename=outpdf,plot=meta_plot,height=4,width=8,dpi=1200)
ggsave(filename=outpng,plot=meta_plot,height=4,width=8,dpi=1200)

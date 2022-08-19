log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log,type="output")
sink(log,type="message")

library("Rcpp")
library("tidyverse")
library("reshape2")
library("tools")

dir <- as.character(snakemake@params[["dir"]])
dir.create(dir)

# Import scaled matrices and sample names as lists and samples config table
sample_table <- read.table(snakemake@input[["sample_table"]],header=T,sep='\t')
matrices_list <- snakemake@input[["matrices"]]
matrices_read <- lapply(matrices_list,function(x) {read.table(x,sep='\t',header=FALSE,skip=1,row.names=3,check.names=FALSE)})
matrices_df <- lapply(matrices_read, function(x) {x <- x[3:ncol(x)]})
matrices_df <- lapply(matrices_df, function(x) {x <- x[rowSums( x[,-1] ) >= 0,]})
sample_names <- lapply(matrices_list, function(x) {readLines(x,n=1)})

# Set variables for output
outpdf <- as.character(snakemake@output[["meta_profile_pdf"]])
outpng <- as.character(snakemake@output[["meta_profile_png"]])

# Import matrix parameters
feature <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["feature"]])))
score <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["score"]])))
prefix <- gsub("([[:lower:]])([[:upper:]])",perl=TRUE,"\\1 \\2",as.character(snakemake@wildcards[["prefix"]]))
normaliser <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["normaliser"]])))
counting <- gsub("([[:lower:]])([[:upper:]])",perl=TRUE,"\\1 \\2",as.character(snakemake@wildcards[["counts"]])) 
paired <- snakemake@params[["paired"]]
paired

bin_size <- as.numeric(snakemake@params[["bin_size"]])
before <- as.numeric(snakemake@params[["before"]])
start <- as.numeric(snakemake@params[["start"]])
body_length <- as.numeric(snakemake@params[["body_length"]])
end <- as.numeric(snakemake@params[["end"]])
after <- as.numeric(snakemake@params[["after"]])

# Normalise to trimmed size factor, calculate trimmed_means of each coordinate on each matrix and combine into a melted long table of metadata
trimmed_colmean <- lapply(matrices_df, function(x) {apply(x,2,function(y) {mean(y,na.rm=F,trim=snakemake@params[["trim"]])})})
trimmed_size <- lapply(trimmed_colmean,function(x) {mean(x)})
size_mean <- mean(unlist(trimmed_size))
for (i in 1:length(trimmed_size)) {trimmed_size[[i]] <- trimmed_size[[i]]/size_mean}
names(trimmed_size) <- sample_names

trimmed_mean <- NULL
for(i in 1:length(trimmed_size)){trimmed_mean[[i]] <- trimmed_colmean[[i]]/trimmed_size[[i]]}
trimmed_mean <- data.frame(trimmed_colmean)
colnames(trimmed_mean) <- sample_names
trimmed_mean$Position <- c((1:nrow(trimmed_mean)*bin_size-(before+bin_size)))
meta_mean <- melt(trimmed_mean,id.vars="Position")
colnames(meta_mean) <- c("Position","Sample","value")
meta_mean$Condition <- sample_table$condition[match(meta_mean$Sample,sample_table$sample_name)]
for (i in (1:nrow(meta_mean))) {meta_mean$rep[i] <- paste("Replicate ",which(sample_table$sample_name[sample_table$condition==meta_mean$Condition[i]]==meta_mean$Sample[i]),sep="")} 

head(meta_mean,10)
write.table(data.frame(trimmed_size),snakemake@output[["size_table"]],quote=F,row.names=F,sep='\t',col.names=names(trimmed_size))

# Set meta-profile colours based on sample config table
meta_mean$colour <- sample_table$colour[match(meta_mean$Sample,sample_table$sample_name)]
sample_colours <- as.character(unique(meta_mean$colour))
sample_names <- as.character(unique(meta_mean$Sample))
names(sample_colours) <- sample_names

heat_colours <- c("dodgerblue","blue","navy","black","red","orange","yellow")

# Generate profile plot parameters
max_val <- ifelse(max(meta_mean$value) >= 0, signif(1.1*max(meta_mean$value),2), 0)
min_val <- ifelse(min(meta_mean$value) >= 0, 0 , signif(1.1*min(meta_mean$value),2))

xlim <- c((0-before),(start+body_length+end+after))
xbrks <- c(0,body_length+start+end,if(before>0)(signif(0-(before*1.5),1)/2) else NULL,if(after>0)(body_length+start+end+signif(after*1.5,1)/2) else NULL,if(start>0)start else NULL,if(end>0)(body_length+start) else NULL)

names(xbrks) <- c("TSS","TTS",if(before>0)(paste("-",signif(before*1.5,1)/2,sep="")) else NULL,if(after>0)(paste("+",signif(after*1.5,1)/2,sep="")) else NULL,if(start>0)(paste("+",start,sep="")) else NULL,if(end>0)(paste("-",end,sep="")) else NULL)

gene_number <- paste("n=",format(nrow(matrices_df[[1]]),big.mark=",",scientific=F),sep="")

title <- paste(prefix, " ", feature," Gene Body",sep="")
subtitle <- paste(score, " Normalised to ",normaliser," ",counting," (",gene_number,")",sep="")
# Plot meta-profile

meta_plot <- ggplot(meta_mean,mapping=aes(x=Position,y=value,colour=Sample)) +
geom_line(size=0.5) +
scale_colour_manual(breaks=sample_names,values=sample_colours) +
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
geom_vline(xintercept=(body_length+start+end),color="black",alpha=0.2,linetype=5) +
geom_vline(xintercept=start,color="black",alpha=0.2,linetype=5) +
geom_vline(xintercept=(body_length+start),color="black",alpha=0.2,linetype=5) +
scale_y_continuous(limits=c(min_val,max_val),breaks=c(min_val,min_val/2,0,max_val/2,max_val)) +
geom_hline(yintercept=0,alpha=0.6) +
scale_x_continuous(limits=xlim,breaks=xbrks) +
ggtitle(title, subtitle=subtitle) +
xlab("Relative Position (bps)") +
ylab("Normalised Counts") +
theme(panel.background=element_rect(fill="White",colour="white"), panel.border=element_rect(fill=NA,colour="black",size=0.7), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.7),axis.line.x.top=element_line(colour="black",size=0.7),axis.line.y.right=element_line(colour="black",size=0.7))
ggsave(filename=outpdf,plot=meta_plot,height=4,width=8,dpi=600)
ggsave(filename=outpng,plot=meta_plot,height=4,width=8,dpi=600)

# Plot per experimental condition meta profiles against control
control <- as.character(snakemake@params[["control"]])
experiments <- unique(as.character(sample_table$condition[sample_table$condition != control]))
experiments

for (i in experiments) {

contrast <- paste(toTitleCase(i),"vs",toTitleCase(control))
exp_meta_mean <- meta_mean[meta_mean$Condition==control|meta_mean$Condition==i,]
exp_meta_mean$colour <- sample_table$colour[match(exp_meta_mean$Sample,sample_table$sample_name)]
exp_sample_colours <- as.character(unique(exp_meta_mean$colour))
exp_condition_colours <- as.character(sample_table$colour[match(unique(exp_meta_mean$Condition),sample_table$condition)])
exp_sample_names <- as.character(unique(exp_meta_mean$Sample))
exp_conditions <- as.character(unique(exp_meta_mean$Condition))
names(exp_sample_colours) <- exp_sample_names
names(exp_condition_colours) <- exp_conditions

if (paired) {

exp_meta_plot <- ggplot(exp_meta_mean,mapping=aes(x=Position,y=value,colour=Condition)) +
facet_wrap(vars(rep),ncol=1,strip.position="top") + 
geom_line(size=0.5,stat="summary_bin",binwidth=bin_size) +
scale_colour_manual(breaks=exp_conditions,values=exp_condition_colours) +
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
geom_vline(xintercept=(body_length+start+end),color="black",alpha=0.2,linetype=5) +
geom_vline(xintercept=start,color="black",alpha=0.2,linetype=5) +
geom_vline(xintercept=(body_length+start),color="black",alpha=0.2,linetype=5) +
scale_y_continuous(limits=c(min_val,max_val),breaks=c(min_val,min_val/2,0,max_val/2,max_val)) +
geom_hline(yintercept=0,alpha=0.6) +
scale_x_continuous(limits=xlim,breaks=xbrks) +
ggtitle(paste(contrast, title, sep=" "), subtitle=subtitle) +
xlab("Relative Position (bps)") +
ylab("Normalised Counts") +
theme(panel.background=element_rect(fill="White",colour="white"), strip.text=element_text(face="bold"),strip.background=element_rect(colour="white",fill="white",size=0.1), panel.border=element_rect(fill=NA,colour="black",size=0.7), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1),axis.line.x.top=element_line(colour="black",size=0.1),axis.line.y.right=element_line(colour="black",size=0.1))
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,"_",subtitle,".pdf",sep="")),path=dir,plot=exp_meta_plot,height=(length(experiments))*3,width=8,dpi=600)
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,"_",subtitle,".png",sep="")),path=dir,plot=exp_meta_plot,height=(length(experiments))*3,width=8,dpi=600)

} else { 

exp_meta_plot <- ggplot(exp_meta_mean,mapping=aes(x=Position,y=value,colour=Sample)) + 
geom_line(size=0.5,stat="summary_bin",binwidth=bin_size) + 
scale_colour_manual(breaks=exp_sample_names,values=exp_sample_colours) + 
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) + 
geom_vline(xintercept=(body_length+start+end),color="black",alpha=0.2,linetype=5) + 
geom_vline(xintercept=start,color="black",alpha=0.2,linetype=5) +
geom_vline(xintercept=(body_length+start),color="black",alpha=0.2,linetype=5) +
scale_y_continuous(limits=c(min_val,max_val),breaks=c(min_val,min_val/2,0,max_val/2,max_val)) + 
geom_hline(yintercept=0,alpha=0.6) + 
scale_x_continuous(limits=xlim,breaks=xbrks) + 
ggtitle(paste(contrast, title, sep=" "), subtitle=subtitle) + 
xlab("Relative Position (bps)") + 
ylab("Normalised Counts") + 
theme(panel.background=element_rect(fill="White",colour="white"), panel.border=element_rect(fill=NA,colour="black",size=0.7), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1),axis.line.x.top=element_line(colour="black",size=0.1),axis.line.y.right=element_line(colour="black",size=0.1))
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,"_",subtitle,".pdf",sep="")),path=dir,plot=exp_meta_plot,height=4,width=8,dpi=600)
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,"_",subtitle,".png",sep="")),path=dir,plot=exp_meta_plot,height=4,width=8,dpi=600)
}
}

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log,type="output")
sink(log,type="message")

library("Rcpp")
library("tidyverse")
library("reshape2")
library("tools")

# SETUP
# Set up multi-threading with BiocParallel
parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dir <- as.character(snakemake@params[["dir"]])
dir.create(dir)

# Import matrices and samples config table
sample_table <- read.table(snakemake@input[["sample_table"]],header=T,sep='\t')

sense_matrices_list <- snakemake@input[["sense_matrices"]]
sense_matrices_read <- lapply(sense_matrices_list,function(x) {read.table(x,sep='\t',header=F,skip=1,row.names=3,check.names=FALSE)})
sense_matrices_df <- lapply(sense_matrices_read, function(x) {x[3:ncol(x)]})
sense_matrices_df <- lapply(sense_matrices_df, function(x) {x <- x[rowSums( x[,-1] ) >= 0,]})
sense_names <- lapply(sense_matrices_list, function(x) {readLines(x,n=1)})
names(sense_matrices_df) <- gsub("_sense","",sense_names)

antisense_matrices_list <- snakemake@input[["antisense_matrices"]]
antisense_matrices_read <- lapply(antisense_matrices_list,function(x) {read.table(x,sep='\t',header=F,skip=1,row.names=3,check.names=FALSE)})
antisense_matrices_df <- lapply(antisense_matrices_read, function(x) {x[3:ncol(x)]})
antisense_matrices_df <- lapply(antisense_matrices_df, function(x) {x <- x[rowSums( x[,-1] ) >= 0,]})
antisense_names <- lapply(antisense_matrices_list, function(x) {readLines(x,n=1)})
names(antisense_matrices_df) <- gsub("_antisense","",antisense_names)
lapply(sense_matrices_df,function(x) {head(x[1:15],10)})


# Set variables for output
outpdf <- as.character(snakemake@output[["meta_profile_pdf"]])
outpng <- as.character(snakemake@output[["meta_profile_png"]])

# Import matrix parameters as variables
feature <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["feature"]])))
score <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["score"]])))
prefix <- gsub("([[:lower:]])([[:upper:]])",perl=TRUE,"\\1 \\2",as.character(snakemake@wildcards[["prefix"]]))
normaliser <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["normaliser"]])))
counting <- gsub("([[:lower:]])([[:upper:]])",perl=TRUE,"\\1 \\2",as.character(snakemake@wildcards[["counts"]]))
paired <- snakemake@params[["paired"]]
paired

bin_size <- as.numeric(snakemake@params[["bin_size"]])
prompt <- as.numeric(snakemake@params[["prompt"]])
tss <- as.numeric(snakemake@params[["tss"]])

# META-PROFILES
# Normalise to sense strand total size, calculate trimmed means of each coordinate on each matrix and combine into a melted long table of metadata
sense_trimmed_colmean <- lapply(sense_matrices_df, function(x) {apply(x,2,function(y) {mean(y,na.rm=F,trim=snakemake@params[["trim"]])})})
size_factors <- read.table(snakemake@input[["size_table"]],header=T,sep='\t',check.names=FALSE)
lapply(sense_trimmed_colmean,function(x) {head(x,10)})
size_factors

sense_trimmed_mean <- NULL
for(i in 1:length(sense_trimmed_colmean)){sense_trimmed_mean[[i]] <- sense_trimmed_colmean[[i]]/size_factors[1,(names(size_factors) == names(sense_trimmed_colmean)[[i]])]}
sense_trimmed_mean <- data.frame(sense_trimmed_mean)
colnames(sense_trimmed_mean) <- sense_names
sense_trimmed_mean$Position <- c((1:nrow(sense_trimmed_mean)*bin_size-(prompt+bin_size)))
sense_meta_mean <- melt(sense_trimmed_mean,id.vars="Position")
sense_meta_mean$Sample <- gsub("_sense","",sense_meta_mean$variable)

antisense_trimmed_colmean <- lapply(antisense_matrices_df, function(x) {apply(x,2,function(y) {mean(y,na.rm=F,trim=snakemake@params[["trim"]])})})
antisense_trimmed_mean <- NULL
for(i in 1:length(antisense_trimmed_colmean)){antisense_trimmed_mean[[i]] <- antisense_trimmed_colmean[[i]]/size_factors[1,(names(size_factors) == names(antisense_trimmed_colmean)[[i]])]}
antisense_trimmed_mean <- data.frame(antisense_trimmed_colmean)
colnames(antisense_trimmed_mean) <- antisense_names
antisense_trimmed_mean$Position <- c(1:nrow(antisense_trimmed_mean)*bin_size-(prompt+bin_size))
antisense_meta_mean <- melt(antisense_trimmed_mean,id.vars="Position")
antisense_meta_mean$value <- 0-antisense_meta_mean$value
antisense_meta_mean$Sample <- gsub("_antisense","",antisense_meta_mean$variable)

meta_mean <- rbind(sense_meta_mean,antisense_meta_mean)
meta_mean$Condition <- sample_table$condition[match(meta_mean$Sample,sample_table$sample_name)]
for (i in (1:nrow(meta_mean))) {meta_mean$rep[i] <- paste("Replicate ",which(sample_table$sample_name[sample_table$condition==meta_mean$Condition[i]]==meta_mean$Sample[i]),sep="")}

head(meta_mean,10)


# Set meta-profile colours based on sample config table & heatmap colour theme
meta_mean$colour <- sample_table$colour[match(meta_mean$Sample,sample_table$sample_name)]
sample_colours <- as.character(unique(meta_mean$colour))
sample_names <- as.character(unique(meta_mean$Sample))
sample_colours
sample_names
names(sample_colours) <- sample_names
sample_colours

heat_colours <- c("dodgerblue","blue","navy","black","red","orange","yellow")

# Generate profile plot parameters
max_val <- ifelse(max(meta_mean$value) >= 0, signif(1.1*max(meta_mean$value),2), 0)
min_val <- ifelse(min(meta_mean$value) >= 0, 0 , signif(1.1*min(meta_mean$value),2))
xlim <- c((0-prompt),(tss))
xbrks <- c(0,if(prompt>0)(signif(0-(prompt*1.5),1)/2) else NULL,if(tss>0)(signif(tss*1.5,1)/2) else NULL, if(tss>0)(signif(tss*1.5,1)/4) else NULL)
names(xbrks) <- c("TSS",if(prompt>0)(paste("-",signif(prompt*1.5,1)/2,sep="")) else NULL,if(tss>0)(paste("+",signif(tss*1.5,1)/2,sep="")) else NULL,if(tss>0)(paste("+",signif(tss*1.5,1)/4,sep="")) else NULL)
gene_number <- paste("n=",format(nrow(sense_matrices_df[[1]]),big.mark=",",scientific=F),sep="")

title <- paste(prefix, " ", feature," PROMPT-TSS",sep="")
subtitle <- paste(score, " Normalised to ",normaliser," ",counting," (",gene_number,")",sep="")

# Plot meta-profile per sample
meta_plot <- ggplot(meta_mean,mapping=aes(x=Position,y=value,group=variable,colour=Sample)) + 
geom_line(size=0.5,stat="summary_bin",binwidth=bin_size) + 
scale_colour_manual(breaks=sample_names,values=sample_colours) + 
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) + 
scale_y_continuous(limits=c(min_val,max_val),breaks=c(min_val,min_val/2,0,max_val/2,max_val)) + 
geom_hline(yintercept=0,alpha=0.6) + 
scale_x_continuous(limits=xlim,breaks=xbrks) + 
ggtitle(title, subtitle=subtitle) +
xlab("Relative Position (bps)") + 
ylab("Normalised Counts") + 
theme(panel.background=element_rect(fill="White",colour="white"),panel.border=element_rect(fill=NA,colour="black",size=0.7), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.7),axis.line.x.top=element_line(colour="black",size=0.7),axis.line.y.right=element_line(colour="black",size=0.7))
ggsave(filename=outpdf,plot=meta_plot,height=4,width=8,dpi=1200)
ggsave(filename=outpng,plot=meta_plot,height=4,width=8,dpi=1200)

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

exp_meta_plot <- ggplot(exp_meta_mean,mapping=aes(x=Position,y=value,group=variable,colour=Condition)) +
facet_wrap(vars(rep),ncol=1,strip.position="top") +
geom_line(size=0.5,stat="summary_bin",binwidth=bin_size) +
scale_colour_manual(breaks=exp_conditions,values=exp_condition_colours) +
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
scale_y_continuous(limits=c(min_val,max_val),breaks=c(min_val,min_val/2,0,max_val/2,max_val)) +
geom_hline(yintercept=0,alpha=0.6) +
scale_x_continuous(limits=xlim,breaks=xbrks) +
ggtitle(paste(contrast, title, sep=" "), subtitle=subtitle) +
xlab("Relative Position (bps)") +
ylab("Normalised Counts") +
theme(panel.background=element_rect(fill="White",colour="white"), strip.text=element_text(face="bold"),strip.background=element_rect(colour="white",fill="white",size=0.1), panel.border=element_rect(fill=NA,colour="black",size=0.7), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1),axis.line.x.top=element_line(colour="black",size=0.1),axis.line.y.right=element_line(colour="black",size=0.1))
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,".pdf",sep="")),path=dir,plot=exp_meta_plot,height=(length(experiments))*3,width=8,dpi=600)
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,".png",sep="")),path=dir,plot=exp_meta_plot,height=(length(experiments))*3,width=8,dpi=600)

} else {

exp_meta_plot <- ggplot(exp_meta_mean,mapping=aes(x=Position,y=value,group=variable,colour=Sample)) +
geom_line(size=0.5,stat="summary_bin",binwidth=bin_size) +
scale_colour_manual(breaks=exp_sample_names,values=exp_sample_colours) +
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
scale_y_continuous(limits=c(min_val,max_val),breaks=c(min_val,min_val/2,0,max_val/2,max_val)) +
geom_hline(yintercept=0,alpha=0.6) +
scale_x_continuous(limits=xlim,breaks=xbrks) +
ggtitle(paste(contrast, title, sep=" "), subtitle=subtitle) +
xlab("Relative Position (bps)") +
ylab("Normalised Counts") +
theme(panel.background=element_rect(fill="White",colour="white"), panel.border=element_rect(fill=NA,colour="black",size=0.7), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1),axis.line.x.top=element_line(colour="black",size=0.1),axis.line.y.right=element_line(colour="black",size=0.1))
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,".pdf",sep="")),path=dir,plot=exp_meta_plot,height=4,width=8,dpi=600)
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,".png",sep="")),path=dir,plot=exp_meta_plot,height=4,width=8,dpi=600)
}
}


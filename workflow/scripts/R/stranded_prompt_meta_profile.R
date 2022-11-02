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
biotype <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["biotype"]])))
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
sense_trimmed_mean <- data.frame(sense_trimmed_colmean)
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

heat_colours <- c("dodgerblue","blue","black","red","orange")

# Generate profile plot parameters
max_val <- ifelse(max(meta_mean$value) >= 0, signif(1.1*max(meta_mean$value),2), 0)
min_val <- ifelse(min(meta_mean$value) >= 0, 0 , signif(1.1*min(meta_mean$value),2))
xlim <- c((0-prompt),(tss))
xbrks <- c(0,if(prompt>0)(signif(0-(prompt*1.5),1)/2) else NULL,if(tss>0)(signif(tss*1.5,1)/2) else NULL, if(tss>0)(signif(tss*1.5,1)/4) else NULL)
names(xbrks) <- c("TSS",if(prompt>0)(paste("-",signif(prompt*1.5,1)/2,sep="")) else NULL,if(tss>0)(paste("+",signif(tss*1.5,1)/2,sep="")) else NULL,if(tss>0)(paste("+",signif(tss*1.5,1)/4,sep="")) else NULL)
gene_number <- paste("n=",format(nrow(sense_matrices_df[[1]]),big.mark=",",scientific=F),sep="")

title <- paste(prefix, " ", biotype," PROMPT-TSS",sep="")
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

sense_df_cond <- NULL
#for (i in 1:(length(sense_matrices_df))) {sense_df_cond[[i]] <- sense_matrices_df[[i]]/size_factors[1,(names(size_factors) == names(sense_matrices_df)[[i]])]}
for (i in 1:(length(sense_matrices_df))) {sense_df_cond[[i]] <- sense_matrices_df[[i]]}
names(sense_df_cond) <- sample_table$condition[match(names(sense_matrices_df),sample_table$sample_name)]
for (i in 1:(length(sense_df_cond))) {sense_df_cond[[i]]$Rank <- nrow(sense_df_cond[[i]]):1}

antisense_df_cond <- NULL
#for (i in 1:(length(antisense_matrices_df))) {antisense_df_cond[[i]] <- antisense_matrices_df[[i]]/size_factors[1,(names(size_factors) == names(sense_matrices_df)[[i]])]}
for (i in 1:(length(antisense_matrices_df))) {antisense_df_cond[[i]] <- antisense_matrices_df[[i]]}
names(antisense_df_cond) <- sample_table$condition[match(names(antisense_matrices_df),sample_table$sample_name)]
for (i in 1:(length(antisense_df_cond))) {antisense_df_cond[[i]]$Rank <- nrow(antisense_df_cond[[i]]):1}

sense_control_df <- sense_matrices_df[names(sense_df_cond) == control]
sense_control_mean <- NULL
for (a in (1:length(sense_control_df))) {sense_control_mean[[a]] <- sense_control_df[[a]]}
for (a in (1:length(sense_control_mean))) {sense_control_mean[[a]]$Rank <- nrow(sense_control_mean[[a]]):1}
sense_control_mean <- do.call(rbind,sense_control_mean)
sense_control_mean <- aggregate(sense_control_mean[,1:(ncol(sense_control_mean)-1)],by=list(sense_control_mean$Rank),mean)
sense_control_mean <- sense_control_mean[2:(ncol(sense_control_mean)-1)]
mean_add <- 1/mean(unlist(size_factors[1,]))

antisense_control_df <- antisense_matrices_df[names(antisense_df_cond) == control]
antisense_control_mean <- NULL
for (a in (1:length(antisense_control_df))) {antisense_control_mean[[a]] <- antisense_control_df[[a]]}
for (a in (1:length(antisense_control_mean))) {antisense_control_mean[[a]]$Rank <- nrow(antisense_control_mean[[a]]):1}
antisense_control_mean <- do.call(rbind,antisense_control_mean)
antisense_control_mean <- aggregate(antisense_control_mean[,1:(ncol(antisense_control_mean)-1)],by=list(antisense_control_mean$Rank),mean)
antisense_control_mean <- antisense_control_mean[2:(ncol(antisense_control_mean)-1)]

file_name <- paste(biotype,"_promptTSS_normalised_",score,"_heatmap",sep="")

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
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,"_",subtitle,".pdf",sep="")),path=dir,plot=exp_meta_plot,height=(length(unique(exp_meta_mean$rep)))*3,width=8,dpi=600)
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,"_",subtitle,".png",sep="")),path=dir,plot=exp_meta_plot,height=(length(unique(exp_meta_mean$rep)))*3,width=8,dpi=600)

sense_exp_df <- sense_matrices_df[names(sense_df_cond) == i]
sense_start <- floor((prompt-100)/bin_size)
sense_fold_change <- NULL
for (a in (1:length(sense_exp_df))) {sense_fold_change[[a]] <- log2((sense_exp_df[[a]][sense_start:(ncol(sense_exp_df[[a]]))]+mean_add)/(sense_control_df[[a]][sense_start:(ncol(sense_exp_df[[a]]))]+mean_add))}
names(sense_fold_change) <- names(sense_exp_df)
for (a in (1:length(sense_fold_change))) {sense_fold_change[[a]]$Rank <- rowSums(sense_fold_change[[a]])}
for (a in (1:length(sense_fold_change))) {sense_fold_change[[a]] <- sense_fold_change[[a]][sense_fold_change[[a]]$Rank != 0,]}
for (a in (1:length(sense_fold_change))) {sense_fold_change[[a]]$Rank <- rank(sense_fold_change[[a]]$Rank,ties.method="last")}
#for (a in (1:length(sense_fold_change))) {sense_fold_change[[a]]$Rank <- nrow(sense_fold_change[[a]]):1}
for (a in (1:length(sense_fold_change))) {sense_fold_change[[a]]$rep <- NA}
for (a in (1:length(sense_fold_change))) {sense_fold_change[[a]]$rep <- paste("Replicate", which(sample_table$sample_name[sample_table$condition==i]==names(sense_fold_change)[a]),sep=" ")}
sense_fold_change <- do.call(rbind,sense_fold_change)
colnames(sense_fold_change)[1:(ncol(sense_fold_change)-2)] <-((1:(ncol(sense_fold_change)-2)*bin_size)-floor(100/bin_size)*bin_size)
head(sense_fold_change[,(ncol(sense_fold_change)-10):ncol(sense_fold_change)],10)
sense_fold_change <- melt(sense_fold_change,id.vars=c("Rank","rep"))
colnames(sense_fold_change) <- c("Rank","rep","Position","log2_FC")
sense_fold_change$Position <- as.numeric(sense_fold_change$Position)

sense_xlim <- c(min(sense_fold_change$Position),max(sense_fold_change$Position))
sense_fold_max <- max(max(sense_fold_change$log2_FC),abs(min(sense_fold_change$log2_FC)))
sense_fold_lim <- c(0-sense_fold_max,sense_fold_max)

antisense_exp_df <- antisense_matrices_df[names(antisense_df_cond) == i]
antisense_end <- floor((prompt+100)/bin_size)
antisense_fold_change <- NULL
for (a in (1:length(antisense_exp_df))) {antisense_fold_change[[a]] <- log2((antisense_exp_df[[a]][1:antisense_end]+mean_add)/(antisense_control_df[[a]][1:antisense_end]+mean_add))}
names(antisense_fold_change) <- names(antisense_exp_df)
for (a in (1:length(antisense_fold_change))) {antisense_fold_change[[a]]$Rank <- rowSums(antisense_fold_change[[a]])}
for (a in (1:length(antisense_fold_change))) {antisense_fold_change[[a]] <- antisense_fold_change[[a]][antisense_fold_change[[a]]$Rank != 0,]}
for (a in (1:length(antisense_fold_change))) {antisense_fold_change[[a]]$Rank <- rank(antisense_fold_change[[a]]$Rank,ties.method="last")}
#for (a in (1:length(antisense_fold_change))) {antisense_fold_change[[a]]$Rank <- nrow(antisense_fold_change[[a]]):1}
for (a in (1:length(antisense_fold_change))) {antisense_fold_change[[a]]$rep <- NA}
for (a in (1:length(antisense_fold_change))) {antisense_fold_change[[a]]$rep <- paste("Replicate", which(sample_table$sample_name[sample_table$condition==i]==names(antisense_fold_change)[a]),sep=" ")}
antisense_fold_change <- do.call(rbind,antisense_fold_change)
colnames(antisense_fold_change)[1:(ncol(antisense_fold_change)-2)] <-((1:(ncol(antisense_fold_change)-2)*bin_size)-floor(prompt/bin_size)*bin_size)
head(antisense_fold_change[,(ncol(antisense_fold_change)-10):ncol(antisense_fold_change)],10)
antisense_fold_change <- melt(antisense_fold_change,id.vars=c("Rank","rep"))
colnames(antisense_fold_change) <- c("Rank","rep","Position","log2_FC")
antisense_fold_change$Position <- as.numeric(antisense_fold_change$Position)

antisense_xlim <- c(min(antisense_fold_change$Position),max(antisense_fold_change$Position))
antisense_fold_max <- max(max(antisense_fold_change$log2_FC),abs(min(antisense_fold_change$log2_FC)))
antisense_fold_lim <- c(0-antisense_fold_max,antisense_fold_max)

# Plot heatmaps
sense_heatmap <- ggplot(sense_fold_change,aes(x=Position,y=Rank,fill=log2_FC)) +
facet_wrap(vars(rep),nrow=1,strip.position="top") +
geom_raster() +
scale_fill_gradientn(colours = heat_colours,limits=sense_fold_lim) +
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
scale_y_discrete(limits=as.integer(c(1,max(sense_fold_change$Rank))),breaks=NULL) +
scale_x_discrete(limits=as.integer(sense_xlim)) +
ggtitle(paste(as.character(i)," vs ",as.character(control)," Sense ",title,sep="")) +
xlab("Relative Position (bps)") +
ylab("Total log2 Fold Change") +
theme(strip.text=element_text(face="bold"),strip.background=element_rect(colour="white",fill="white",size=0.1), panel.background=element_rect(fill="White",colour="white"),panel.border=element_rect(fill=NA,colour="black",size=0.7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1))

antisense_heatmap <- ggplot(antisense_fold_change,aes(x=Position,y=Rank,fill=log2_FC)) +
facet_wrap(vars(rep),nrow=1,strip.position="top") +
geom_raster() +
scale_fill_gradientn(colours = heat_colours,limits=antisense_fold_lim) +
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
scale_y_discrete(limits=as.integer(c(1,ncol(sense_fold_change)-1)),breaks=NULL) +
scale_x_discrete(limits=as.integer(antisense_xlim)) +
ggtitle(paste(as.character(i)," vs ",as.character(control)," Antisense ",title,sep="")) +
xlab("Relative Position (bps)") +
ylab("Total log2 Fold Change") +
theme(strip.text=element_text(face="bold"),strip.background=element_rect(colour="white",fill="white",size=0.1), panel.background=element_rect(fill="White",colour="white"), panel.border=element_rect(fill=NA,colour="black",size=0.7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1))


ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_sense_",file_name,".pdf",sep=""),path=dir,plot=sense_heatmap,height=11,width=6,dpi=600)
ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_sense_",file_name,".png",sep=""),path=dir,plot=sense_heatmap,height=11,width=6,dpi=600)

ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_antisense_",file_name,".pdf",sep=""),path=dir,plot=antisense_heatmap,height=11,width=6,dpi=600)
ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_antisense_",file_name,".png",sep=""),path=dir,plot=antisense_heatmap,height=11,width=6,dpi=600)

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
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,"_",subtitle,".pdf",sep="")),path=dir,plot=exp_meta_plot,height=4,width=8,dpi=600)
ggsave(filename=gsub(" ","_",paste(contrast,"_",title,"_",subtitle,".png",sep="")),path=dir,plot=exp_meta_plot,height=4,width=8,dpi=600)

sense_exp_df <- sense_df_cond[names(sense_df_cond) == i]
sense_exp_df <- do.call(rbind,sense_exp_df)
sense_exp_mean <- aggregate(sense_exp_df,by=list(sense_exp_df$Rank),mean)
sense_exp_mean <- sense_exp_mean[2:(ncol(sense_exp_mean)-1)]
sense_fold_change <-log2(((sense_exp_mean+mean_add)/(sense_control_mean+mean_add)))
colnames(sense_fold_change) <-((0-prompt):(ncol(sense_fold_change)-prompt-1))
sense_fold_change$Rank <- rowSums(sense_fold_change)
sense_fold_change <- sense_fold_change[sense_fold_change$Rank != 0,]
sense_fold_change$Rank <- rank(sense_fold_change$Rank,ties.method="last")
#sense_fold_change$Rank <- nrow(sense_fold_change):1
sense_fold_change <- melt(sense_fold_change,id.vars="Rank")
colnames(sense_fold_change) <- c("Rank","Position","log2_FC")
sense_fold_change$Position <- as.numeric(sense_fold_change$Position)

sense_fold_max <- max(max(sense_fold_change$log2_FC),abs(min(sense_fold_change$log2_FC)))
sense_fold_lim <- c(0-sense_fold_max,sense_fold_max)

antisense_exp_df <- antisense_df_cond[names(antisense_df_cond) == i]
antisense_exp_df <- do.call(rbind,antisense_exp_df)
antisense_exp_mean <- aggregate(antisense_exp_df,by=list(antisense_exp_df$Rank),mean)
antisense_exp_mean <- antisense_exp_mean[2:(ncol(antisense_exp_mean)-1)]
antisense_fold_change <-log2(((antisense_exp_mean+mean_add)/(antisense_control_mean+mean_add)))
colnames(antisense_fold_change) <- ((0-prompt):(ncol(antisense_fold_change)-prompt-1))
antisense_fold_change$Rank <- rowSums(antisense_fold_change)
antisense_fold_change <- antisense_fold_change[antisense_fold_change$Rank != 0,]
antisense_fold_change$Rank <- rank(antisense_fold_change$Rank,ties.method="last")
#antisense_fold_change$Rank <- nrow(antisense_fold_change):1
antisense_fold_change <- melt(antisense_fold_change,id.vars="Rank")
colnames(antisense_fold_change) <- c("Rank","Position","log2_FC")
antisense_fold_change$Position <- as.numeric(antisense_fold_change$Position)

# Plot heatmaps
sense_heatmap <- ggplot(sense_fold_change,aes(x=Position,y=Rank,fill=log2_FC)) +
geom_raster() +
scale_fill_gradientn(colours = heat_colours,limits=sense_fold_lim) +
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
scale_y_discrete(limits=as.integer(c(1,ncol(sense_fold_change)-1)),breaks=NULL) +
scale_x_discrete(limits=as.integer(xlim),breaks=xbrks) +
ggtitle(paste(as.character(i)," vs ",as.character(control)," Sense ",title,sep="")) +
xlab("Relative Position (bps)") +
ylab("Total log2 Fold Change") +
theme(panel.background=element_rect(fill="White",colour="white"),panel.border=element_rect(fill=NA,colour="black",size=0.7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1))

antisense_heatmap <- ggplot(antisense_fold_change,aes(x=Position,y=Rank,fill=log2_FC)) +
geom_raster() +
scale_fill_gradientn(colours = heat_colours,limits=antisense_fold_lim) +
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
scale_y_discrete(limits=as.integer(c(1,ncol(sense_fold_change)-1)),breaks=NULL) +
scale_x_discrete(limits=as.integer(xlim),breaks=xbrks) +
ggtitle(paste(as.character(i)," vs ",as.character(control)," Antisense ",title,sep="")) +
xlab("Relative Position (bps)") +
ylab("Total log2 Fold Change") +
theme(panel.background=element_rect(fill="White",colour="white"), panel.border=element_rect(fill=NA,colour="black",size=0.7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1))


ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_sense_",file_name,".pdf",sep=""),path=dir,plot=sense_heatmap,height=5,width=3,dpi=1200)
ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_sense_",file_name,".png",sep=""),path=dir,plot=sense_heatmap,height=5,width=3,dpi=1200)

ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_antisense_",file_name,".pdf",sep=""),path=dir,plot=antisense_heatmap,height=5,width=3,dpi=1200)
ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_antisense_",file_name,".png",sep=""),path=dir,plot=antisense_heatmap,height=5,width=3,dpi=1200)
}
}


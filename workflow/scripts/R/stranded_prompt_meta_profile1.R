log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log,type="output")
sink(log,type="message")

library("Rcpp")
library("tidyverse")
library("reshape2")
library("tools")

# SETUP
# Import matrices and samples config table
sample_table <- read.table(snakemake@input[["sample_table"]],header=T,sep='\t')

sense_matrices_list <- snakemake@input[["sense_matrices"]]
sense_matrices_df <- lapply(sense_matrices_list,function(x) {read.table(x,sep='\t',header=F,skip=1,row.names=2)})
sense_matrices_df <- lapply(sense_matrices_df, function(x) {x[2:ncol(x)]})
sense_names <- lapply(sense_matrices_list, function(x) {readLines(x,n=1)})
names(sense_matrices_df) <- gsub("_sense","",sense_names)

antisense_matrices_list <- snakemake@input[["antisense_matrices"]]
antisense_matrices_df <- lapply(antisense_matrices_list,function(x) {read.table(x,sep='\t',header=F,skip=1,row.names=2)})
antisense_matrices_df <- lapply(antisense_matrices_df, function(x) {x[2:ncol(x)]})
antisense_names <- lapply(antisense_matrices_list, function(x) {readLines(x,n=1)})
names(antisense_matrices_df) <- gsub("_antisense","",antisense_names)

# Set variables for output
outpdf <- snakemake@output[["meta_profile_pdf"]]
outpng <- snakemake@output[["meta_profile_png"]]

# Import matrix parameters as variables
feature <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["feature"]])))
score <- toTitleCase(gsub("_"," ",as.character(snakemake@wildcards[["score"]])))
bin_size <- as.numeric(snakemake@params[["bin_size"]])
prompt <- as.numeric(snakemake@params[["prompt"]])
tss <- as.numeric(snakemake@params[["tss"]])

# META-PROFILES
# Normalise to sense strand total size, calculate trimmed means of each coordinate on each matrix and combine into a melted long table of metadata
sense_trimmed_colmean <- lapply(sense_matrices_df, function(x) {apply(x,2,function(y) {mean(y,na.rm=F,trim=snakemake@params[["trim"]])})})
#sense_trimmed_size <- lapply(sense_trimmed_colmean, function(x) {mean(x[1:(prompt-200)])})

#sense_trimmed_mean <- NULL
#for(i in 1:length(sense_trimmed_size)){sense_trimmed_mean[[i]] <- sense_trimmed_colmean[[i]]/sense_trimmed_size[[i]]}
sense_trimmed_mean <- data.frame(sense_trimmed_colmean)
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


# Set meta-profile colours based on sample config table & heatmap colour theme
meta_mean$colour <- sample_table$colour[match(meta_mean$Sample,sample_table$sample_name)]
sample_colours <- as.character(unique(meta_mean$colour))
sample_names <- as.character(unique(meta_mean$Sample))
sample_colours
sample_names
names(sample_colours) <- sample_names
sample_colours

heat_colours <- c("blue","navy","black","red","yellow")

# Generate profile plot parameters
max_val <- ifelse(max(meta_mean$value) >= 0, signif(1.1*max(meta_mean$value),2), 0)
min_val <- ifelse(min(meta_mean$value) >= 0, 0 , signif(1.1*min(meta_mean$value),2))
xlim <- c((0-prompt),(tss))
xbrks <- c(0,if(prompt>0)(signif(0-(prompt*1.5),1)/2) else NULL,if(tss>0)(signif(tss*1.5,1)/2) else NULL, if(tss>0)(signif(tss*1.5,1)/4) else NULL)
names(xbrks) <- c("TSS",if(prompt>0)(paste("-",signif(prompt*1.5,1)/2,sep="")) else NULL,if(tss>0)(paste("+",signif(tss*1.5,1)/2,sep="")) else NULL,if(tss>0)(paste("+",signif(tss*1.5,1)/4,sep="")) else NULL)
gene_number <- paste("n=",format(nrow(sense_matrices_df[[1]]),big.mark=",",scientific=F),sep="")
title <- paste(feature," PROMPT-TSS (",gene_number,") Nomalised ",score,sep="")

# Plot meta-profile per sample
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

# HEATMAP
# Define control and experiment conditions, normalise and name matrices by conditions
control <- as.character(snakemake@params[["control"]])
experiments <- unique(as.character(sample_table$condition[sample_table$condition != control]))
control
experiments

sense_df_cond <- NULL
for (i in 1:(length(sense_trimmed_size))) {sense_df_cond[[i]] <- sense_matrices_df[[i]]/sense_trimmed_size[[i]]}
names(sense_df_cond) <- sample_table$condition[match(names(sense_matrices_df),sample_table$sample_name)]
for (i in 1:(length(sense_df_cond))) {sense_df_cond[[i]]$Rank <- nrow(sense_df_cond[[i]]):1}

antisense_df_cond <- NULL
for (i in 1:(length(sense_trimmed_size))) {antisense_df_cond[[i]] <- antisense_matrices_df[[i]]/sense_trimmed_size[[i]]}
names(antisense_df_cond) <- sample_table$condition[match(names(antisense_matrices_df),sample_table$sample_name)]
for (i in 1:(length(antisense_df_cond))) {antisense_df_cond[[i]]$Rank <- nrow(antisense_df_cond[[i]]):1}

# Generate summary meta data for control
sense_control_df <- sense_df_cond[names(sense_df_cond) == control]
sense_control_df <- do.call(rbind,sense_control_df)
sense_control_mean <- aggregate(sense_control_df,by=list(sense_control_df$Rank),mean)
sense_control_mean <- sense_control_mean[2:(ncol(sense_control_mean)-1)]
mean_add <- 1/mean(unlist(sense_trimmed_size))

antisense_control_df <- antisense_df_cond[names(antisense_df_cond) == control]
antisense_control_df <- do.call(rbind,antisense_control_df)
antisense_control_mean <- aggregate(antisense_control_df,by=list(antisense_control_df$Rank),mean)
antisense_control_mean <- antisense_control_mean[2:(ncol(antisense_control_mean)-1)]

file_name <- paste(as.character(snakemake@wildcards[["feature"]]),"_promptTSS_normalised_",as.character(snakemake@wildcards[["score"]]),"_heatmap",sep="")

# For each expeirmental condtion, generate meta data and fold-change matrices
for (i in experiments) {
sense_exp_df <- sense_df_cond[names(sense_df_cond) == i]
sense_exp_df <- do.call(rbind,sense_exp_df)
sense_exp_mean <- aggregate(sense_exp_df,by=list(sense_exp_df$Rank),mean)
sense_exp_mean <- sense_exp_mean[2:(ncol(sense_exp_mean)-1)]
sense_fold_change <-log2(((sense_exp_mean+mean_add)/(sense_control_mean+mean_add)))
colnames(sense_fold_change) <-((0-prompt):(ncol(sense_fold_change)-prompt-1))
sense_fold_change$Rank <- rowSums(sense_fold_change)
sense_fold_change$Rank <- rank(sense_fold_change$Rank)
sense_fold_change <- melt(sense_fold_change,id.vars="Rank")
colnames(sense_fold_change) <- c("Rank","Position","log2_FC")
sense_fold_change$Position <- as.numeric(sense_fold_change$Position)

antisense_exp_df <- antisense_df_cond[names(antisense_df_cond) == i]
antisense_exp_df <- do.call(rbind,antisense_exp_df)
antisense_exp_mean <- aggregate(antisense_exp_df,by=list(antisense_exp_df$Rank),mean)
antisense_exp_mean <- antisense_exp_mean[2:(ncol(antisense_exp_mean)-1)]
antisense_fold_change <-log2(((antisense_exp_mean+mean_add)/(antisense_control_mean+mean_add)))
colnames(antisense_fold_change) <- ((0-prompt):(ncol(antisense_fold_change)-prompt-1))
antisense_fold_change$Rank <- rowSums(antisense_fold_change)
antisense_fold_change$Rank <- rank(antisense_fold_change$Rank)
antisense_fold_change <- melt(antisense_fold_change,id.vars="Rank")
colnames(antisense_fold_change) <- c("Rank","Position","log2_FC")
antisense_fold_change$Position <- as.numeric(antisense_fold_change$Position)

write.table(sense_fold_change,file=paste(as.character(snakemake@output[["heatmap_dir"]]),"/",as.character(i),"_vs_",as.character(control),"_sense_",file_name,".tab",sep=""),sep='\t',row.names=F,quote=F)
write.table(antisense_fold_change,file=paste(as.character(snakemake@output[["heatmap_dir"]]),"/",as.character(i),"_vs_",as.character(control),"_antisense_",file_name,".tab",sep=""),sep='\t',row.names=F,quote=F)

# Plot heatmaps
sense_heatmap <- ggplot(sense_fold_change,aes(x=Position,y=Rank,fill=log2_FC)) + 
geom_tile() + 
scale_fill_gradientn(colours = heat_colours) + 
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
scale_y_discrete(limits=c(1,ncol(sense_fold_change)-1),breaks=NULL) + 
scale_x_discrete(limits=xlim,breaks=xbrks) +
ggtitle(paste(as.character(i)," vs ",as.character(control)," Sense ",title,sep="")) + 
xlab("Relative Position (bps)") +
ylab("Total log2 Fold Change") +
theme(panel.background=element_rect(fill="White",colour="white"),panel.border=element_rect(fill=NA,colour="black",size=0.7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1))

antisense_heatmap <- ggplot(antisense_fold_change,aes(x=Position,y=Rank,fill=log2_FC)) +
geom_tile() +
scale_fill_gradientn(colours = heat_colours) +
geom_vline(xintercept=0,linetype=5,colour="black",alpha=0.2) +
scale_y_discrete(limits=c(1,ncol(sense_fold_change)-1),breaks=NULL) +
scale_x_discrete(limits=xlim,breaks=xbrks) +
ggtitle(paste(as.character(i)," vs ",as.character(control)," Antisense ",title,sep="")) + 
xlab("Relative Position (bps)") +
ylab("Total log2 Fold Change") +
theme(panel.background=element_rect(fill="White",colour="white"), panel.border=element_rect(fill=NA,colour="black",size=0.7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1))


ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_sense_",file_name,".pdf",sep=""),path=as.character(snakemake@output[["heatmap_dir"]]),plot=sense_heatmap,height=12,width=8,dpi=1200)
ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_sense_",file_name,".png",sep=""),path=as.character(snakemake@output[["heatmap_dir"]]),plot=sense_heatmap,height=12,width=8,dpi=1200)

ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_antisense_",file_name,".pdf",sep=""),path=as.character(snakemake@output[["heatmap_dir"]]),plot=antisense_heatmap,height=12,width=8,dpi=1200)
ggsave(filename=paste(as.character(i),"_vs_",as.character(control),"_antisense_",file_name,".png",sep=""),path=as.character(snakemake@output[["heatmap_dir"]]),plot=antisense_heatmap,height=12,width=8,dpi=1200)
}

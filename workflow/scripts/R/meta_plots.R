
# Process matrices to trimmed mean and/or/maybe spline
for ( i in sense_dirs ) {
mx <- get(paste(i,"_mx",sep=""))
mx <- lapply(mx,function(x) {x[rownames(x) %in% expr_i$featureID,]} )
meta_gene_n <- nrow(mx[[1]])
mx_mean <- data.frame(lapply(mx, function(x) {apply(x,2,function(y) {mean(y,na.rm=F,trim=meta_trim)})}),check.names=F)
mx_data <- data.frame(
  unlist(
    lapply(colnames(mx_mean), function (x) {
      ifelse(i=="sense", mx_mean[paste(x)], 0-mx_mean[paste(x)])
    })),
  unlist(
    lapply(colnames(mx_mean), function (x) {
      replicate(nrow(mx_mean),paste(x))
    })),
  unlist(
    lapply(colnames(mx_mean), function (x) {
      ifelse(i=="sense",as.numeric(rownames(mx_mean)), rev(as.numeric(rownames(mx_mean))))
    })),
  unlist(
    lapply(colnames(mx_mean), function (x) {
      replicate(nrow(mx_mean),paste(x,i,sep="_"))
    }))  
)
names(mx_data) <- c("value","Sample","Position","variable")
assign(paste(i,"_mx_data",sep=""),mx)
}

head(sense_mx_data,10)

if (length(antisense_ls)>0) {
sense_mx_data <- rbind(sense_mx_data,antisense_mx_data)
}

sense_mx_data$Condition <- sample_table$condition[match(sense_mx_data$Sample,sample_table$sample_name)]
sense_mx_data$colour <- sample_table$colour[match(sense_mx_data$Sample,sample_table$sample_name)]
sample_colours <- as.character(unique(sense_mx_data$colour))
sample_names <- as.character(unique(sense_mx_data$Sample))
sample_colours
sample_names
names(sample_colours) <- gsub("_"," ",sample_names)
sample_colours

sense_mx_data$rep <- paste("Replicate", sample_table$replicate[match(sense_mx_data$Sample,sample_table$sample_name)])

max_val <- ifelse(max(sense_mx_data$value) >= 0, signif(1.1*max(sense_mx_data$value),2), 0)
min_val <- ifelse(min(sense_mx_data$value) >= 0, 0 , signif(1.1*min(sense_mx_data$value),2))
xlim <- c(min(sense_mx_data$Position),max(sense_mx_data$Position))

if (section=="body") {
xbrks <- c(0,main_bin)
names(xbrks) <- c("Start","End")
} else {
xbrks <- c(ifelse(len_bef_n>0,0-bef_bin,NULL),0,ifelse(len_aft_n>0,main_bin,NULL))
names(xbrks) <- c(ifelse(len_bef_n>0,len_bef,NULL),toTitleCase(paste(base,section)),ifelse(len_aft_n>0,len_aft,NULL))
}
sense_mx_data$Sample <- gsub("_"," ",sense_mx_data$Sample)

head(sense_mx_data,10)


heat_colours <- c("dodgerblue","blue","black","red","orange")

meta <- ggplot(sense_mx_data,mapping=aes(x=Position,y=value,group=variable,colour=Sample)) +
  geom_line(
    size=0.5,
    stat="summary_bin",
    binwidth=1
  ) +
  scale_colour_manual(
    breaks=mx_samples,
    values=sample_colours
  ) +
  geom_vline(
    xintercept=0,
    linetype=5,
    colour="black",
    alpha=0.2
  ) +
  scale_y_continuous(
    limits=c(min_val,max_val),
    breaks=c(min_val,min_val/2,0,max_val/2,max_val)
  ) +
  geom_hline(
    yintercept=0,
    alpha=0.6
  ) +
  scale_x_continuous(
    limits=xlim,
    breaks=xbrks
  ) +
  ylab("Normalised Counts") +
  theme(
    panel.background=element_rect(fill="White",colour="white"), 
    panel.border=element_rect(fill=NA,colour="black",size=0.7), 
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",size=0.1)
    legend.background=element_rect(fill="White"), 
    legend.key=element_rect(colour="white",fill="White"), 
    axis.line=element_line(colour="black",size=0.1),
    axis.line.x.top=element_line(colour="black",size=0.1),
    axis.line.y.right=element_line(colour="black",size=0.1)
  )

meta_plot <- meta + 
  facet_wrap(vars(rep),ncol=2,strip.position="top") +
  xlab(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE)) +
  ylab("Normalised Coverage") +  

if (rep_pair) {

meta_plot <- meta_plot + facet_wrap(vars(rep),ncol=1,strip.position="top")

} 

meta_caption <- paste(
  "Normalised base-coverage over ", feature_i, " in ", experiment "." 
, sep="" )

meta_plot_title <- paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE),"Coverage Profiles.")
meta_plot_caption <- meta_caption 



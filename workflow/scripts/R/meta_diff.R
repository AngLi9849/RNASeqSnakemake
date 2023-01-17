#meta_gene_n <- sum(rownames(sig_bg)[ (sig_bg$sig2bg >= sig & sig_bg$bg2sig >= bg) ] %in% expr_i$featureID[expr_i$baseMean >= config_min_mean & expr_i$baseMean >= min_rpkm])
#mx_data <- mx_df[mx_df$i_group == i,]

#mx_data <- heat_data
meta_gene_n <- length(unique(heat_data$featureID))

rep_table <- sample_table[sample_table$condition==treatment,]
rep_table$replicate_name <-  paste("Replicate",rep_table$replicate)
rep_colours <- rep_table$colour
names(rep_colours) <- rep_table$replicate_name
rep_lines <- rep_table$replicate
names(rep_lines) <- rep_table$replicate_name
rep_brks <- rep_table$replicate_name

max_val <- ifelse(max(mx_data$log2FoldChange) >= 0, signif(1.1*max(mx_data$log2FoldChange),2), 0)
min_val <- ifelse(min(mx_data$log2FoldChange) >= 0, 0 , signif(1.1*min(mx_data$log2FoldChange),2))
xlim <- c(min(mx_data$Position),max(mx_data$Position))

exp_col <- condition_col[treat]

#mx_data$replicate <- paste("Replicate",mx_data$replicate)

mx_data$replicate <- factor(mx_data$replicate, levels=c("All",rep_table$replicate_name))
mx_data$Sense <- factor(mx_data$Sense, levels=c("Sense","Antisense"))


mx_data$replicate_sense <- paste(mx_data$replicate, mx_data$Sense)

mx_data$replicate_sense <- factor(
  mx_data$replicate_sense, 
  levels=c(
    paste(c("All",rep_table$replicate_name),"Sense"),
    paste(c("All",rep_table$replicate_name),"Antisense")
  )
)


for (reps in c("All","Reps") ) {

  if (reps=="All") {
    mx_rep_data <- mx_data[mx_data$replicate == "All",]
  } else {
    mx_rep_data <- mx_data[mx_data$replicate != "All",]
  }


meta_base <- ggplot(mx_rep_data,mapping=aes(x=Position,y=log2FoldChange)) +
  geom_vline(
    xintercept=0,
    linetype=5,
    colour="black",
    alpha=0.2
  ) +
  scale_y_continuous(
    limits=c(min_val,max_val),
    breaks=c(min_val,0,max_val)
  ) +
  geom_hline(
    yintercept=0,
    alpha=0.6
  ) +
  scale_x_continuous(
    limits=xlim,
    breaks=heat_xbrks,
  ) +
  xlab(feature) +
  ylab(paste("\U0394","log2(",difference,")",sep="")) +
  theme(
    panel.background=element_rect(fill="White",colour="white"), 
    panel.border=element_rect(fill=NA,colour="black",linewidth=0.7), 
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",linewidth=0.1),
    legend.background=element_rect(fill="White"), 
    legend.key=element_rect(colour="white",fill="White"), 
    legend.position="none",
    axis.text=element_text(colour="black"),
    axis.text.x = if (length(heat_xbrks)>2) (element_text(angle = 45, vjust = 1, hjust=1,colour="black")) else (element_text(colour="black")),
    axis.line=element_line(colour="black",linewidth=0.1),
    axis.line.x.top=element_line(colour="black",linewidth=0.1),
    axis.line.y.right=element_line(colour="black",linewidth=0.1)
  )

if (section=="body") {

meta_base <- meta_base + 
  geom_vline(
    xintercept=main_bin*heat_bin/meta_bin,
    linetype=5,
    colour="black",
    alpha=0.2
  )
}

if (reps=="All") {

meta <- meta_base + 
  geom_smooth(
    linewidth=0.8,
    method="loess",
    n=300,
    span=0.05,
    colour=exp_col,
    alpha=0,
  ) + 
  geom_ribbon(
    aes(
      ymin=log2FoldChange-lfcSE,
      ymax=log2FoldChange+lfcSE
    ),
    fill=exp_col,
    alpha=0.2
  ) +
  scale_colour_manual(
    values=exp_col
  ) +
  labs(
    subtitle=paste(ifelse(common=="","",paste(common," ",sep="")),treat," vs ",control,"\n",protocol_title,"\n","n=",meta_gene_n,sep="")
  ) +
  scale_fill_manual(
    values=exp_col
  )

  if (length(unique(mx_data$Sense)) > 1) {
    meta <- meta + facet_wrap(vars(Sense),ncol=1,strip.position="right")
  }
}

if (reps!="All") {

meta_plot <- meta_base +
  geom_smooth(
    linewidth=0.8,
    method="loess",
    n=300,
    span=0.05,
    alpha=0,
    aes(group=replicate_sense,colour=replicate,linetype=replicate,fill=replicate),
  ) +
  geom_ribbon(
    aes(
      ymin=log2FoldChange-lfcSE,
      ymax=log2FoldChange+lfcSE,
      group=replicate_sense,fill=replicate
    ),
    alpha=0.2
  ) +
  scale_colour_manual(
    breaks=rep_brks,
    values=rep_colours,
  ) +
  scale_fill_manual(
    breaks=rep_brks,
    values=rep_colours,
  ) +
  labs(
    subtitle=paste(ifelse(common=="","",paste(common," ",sep="")), treat," vs ",control,"\n",protocol_title," ",prefix,"\n","n=",meta_gene_n,sep="")
  ) +
  xlab(
    paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE), "Position (bps)")
  ) +
  scale_x_continuous(
    limits=xlim,
    breaks=heat_xbrks,
  ) +  
  theme(
    legend.position = "bottom",
    legend.direction="vertical",
    legend.just="left",
  )

if (rep_pair & length(unique(mx_data$replicate)) > 1) {

rep_lines <- rep(1,nrow(rep_table))
names(rep_lines) <- rep_table$replicate_name
meta_plot <- meta_plot + 
  scale_linetype_manual(
    breaks = rep_brks,
    values = rep_lines 
  ) 

  if (length(unique(mx_data$Sense)) <= 1) {
    meta_plot <- meta_plot + facet_wrap(vars(replicate),ncol=2,strip.position="top")
    meta_plot_h <- min(c(7.2,2.4*ceiling(length(unique(mx_data$replicate))/2)))
  } else {
    meta_plot <- meta_plot + facet_grid(replicate~Sense)
    meta_plot_h <- min(c(7.3,2.4*ceiling(length(unique(mx_data$replicate)))))
  }

} else {

meta_plot <- meta_plot +
  scale_linetype_manual(
    breaks = rep_brks,
    values = rep_lines
  )
  if (length(unique(mx_data$Sense)) <= 1) {
    meta_plot_h <- 4 + 0.1*length(unique(mx_data$replicate))
  } else {
    meta_plot <- meta_plot + facet_wrap(vars(Sense),ncol=1,strip.position="top")
    meta_plot_h <- 6.6 + 0.1*length(unique(mx_data$replicate)) 
  }
}

}
}
meta_caption <- paste(
  "Normalised change in log2 ", difference, " coverage over ",meta_gene_n," ", feature_i, " in ", experiment, "." 
, sep="" )
meta_plot_title <- paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE),"change in log2", difference,  "coverage.")
meta_plot_caption <- meta_caption 

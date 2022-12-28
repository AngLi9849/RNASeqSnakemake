meta_gene_n <- sum(rownames(sig_bg)[ (sig_bg$sig2bg >= sig & sig_bg$bg2sig >= bg) ] %in% expr_i$featureID[expr_i$baseMean >= config_min_mean & expr_i$baseMean >= min_rpkm])
#mx_data <- mx_df[mx_df$i_group == i,]

mx_data <- heat_data
sample_colours <- sample_table$colour
names(sample_colours) <- gsub("_"," ",sample_table$sample_abbr)
sample_lines <- sample_table$replicate
names(sample_lines) <- gsub("_"," ",sample_table$sample_abbr)
sample_brks <- gsub("_"," ",sample_table$sample_abbr)

max_val <- ifelse(max(mx_data$value) >= 0, signif(1.1*max(mx_data$value),2), 0)
min_val <- ifelse(min(mx_data$value) >= 0, 0 , signif(1.1*min(mx_data$value),2))
xlim <- c(min(mx_data$Position),max(mx_data$Position))

mx_data$Sample <- gsub("_"," ",mx_data$Sample)
mx_data$Sample_sense <- paste(mx_data$Sample, mx_data$sense)

head(mx_data,10)

mx_data$Sample <- factor(mx_data$Sample, levels=sample_table$sample_abbr)
mx_data$Sample_sense <- factor(
  mx_data$Sample_sense, 
  levels=c(
    paste(sample_table$sample_abbr[sample_table$cond_abbr == control],"sense"),
    paste(sample_table$sample_abbr[sample_table$cond_abbr == treat],"sense"),
    paste(sample_table$sample_abbr[sample_table$cond_abbr == control],"antisense"),
    paste(sample_table$sample_abbr[sample_table$cond_abbr == treat],"antisense")
  )
)

meta_base <- ggplot(mx_data,mapping=aes(x=Position,y=value)) +
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
    breaks=meta_xbrks,
  ) +
  labs(
    subtitle=paste("n=",meta_gene_n,sep="")
  ) +
  xlab(feature) +
  ylab(paste(meta_y,"Coverage")) +
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
    xintercept=main_bin,
    linetype=5,
    colour="black",
    alpha=0.2
  )
}

meta <- meta_base + 
  geom_smooth(
    linewidth=0.8,
    method="loess",
    n=300,
    span=0.05,
    alpha=0.2,
    aes(group=cond_group,colour=Condition,fill=Condition),
  ) + 
  scale_colour_manual(
    values=condition_col
  ) +
  scale_fill_manual(
    values=condition_col
  )


meta_plot <- meta_base +
  geom_smooth(
    linewidth=0.8,
    method="loess",
    n=300,
    span=0.05,
    alpha=0.2,
    aes(group=Sample_sense,colour=Sample,linetype=Sample,fill=Sample),
  ) +
  scale_colour_manual(
    breaks=sample_brks,
    values=sample_colours,
  ) +
  scale_fill_manual(
    breaks=sample_brks,
    values=sample_colours,
  ) +
  xlab(
    paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE), "Position (bps)")
  ) +
  scale_x_continuous(
    limits=xlim,
    breaks=meta_xbrks,
  ) +  
  ylab(paste("Normalised",meta_y, "Coverage Depth")) + 
  theme(
    legend.position = "right"
  )

if (rep_pair & length(unique(mx_data$rep)) > 1) {

sample_lines <- rep(1,nrow(sample_table))
names(sample_lines) <- gsub("_"," ",sample_table$sample_abbr)
meta_plot <- meta_plot + 
  facet_wrap(vars(rep),ncol=1,strip.position="top") +
  scale_linetype_manual(
    breaks = sample_brks,
    values = sample_lines 
  )

meta_plot_h <- length(unique(mx_data$rep))*min(c(2.3,7/length(unique(mx_data$rep))))

} else {

meta_plot <- meta_plot +
  scale_linetype_manual(
    breaks = sample_brks,
    values = sample_lines
  )
meta_plot_h <- 3.5
  
}

meta_caption <- paste(
  "Normalised coverage over ",meta_gene_n," ", feature_i, " in ", experiment, "." 
, sep="" )
meta_plot_title <- paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE),"Coverage Profiles.")
meta_plot_caption <- meta_caption 



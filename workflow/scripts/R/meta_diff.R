#meta_gene_n <- sum(rownames(sig_bg)[ (sig_bg$sig2bg >= sig & sig_bg$bg2sig >= bg) ] %in% expr_i$featureID[expr_i$baseMean >= config_min_mean & expr_i$baseMean >= min_rpkm])
#mx_data <- mx_df[mx_df$i_group == i,]

mx_data <- heat_data
rep_table <- sample_table[sample_table$condition==treatment,]
rep_table$replicate_name <-  paste("Replicate",rep_table$replicate)
rep_colours <- rep_table$colour
names(rep_colours) <- rep_table$replicate_name
rep_lines <- rep_table$replicate
names(rep_lines) <- rep_table$replicate_name
rep_brks <- rep_table$replicate_name

max_val <- ifelse(max(heat_data$value) >= 0, signif(1.1*max(heat_data$value),2), 0)
min_val <- ifelse(min(mx_data$value) >= 0, 0 , signif(1.1*min(mx_data$value),2))
xlim <- c(min(mx_data$Position),max(mx_data$Position))

exp_col <- condition_col[treat]

#mx_data$replicate <- paste("Replicate",mx_data$replicate)
mx_data$replicate_sense <- paste(mx_data$replicate, mx_data$Sense)

head(mx_data,10)

mx_data$replicate <- factor(mx_data$replicate, levels=c("All",rep_table$replicate_name)
mx_data$replicate_sense <- factor(
  mx_data$replicate_sense, 
  levels=c(
    paste(c("All",rep_table$replicate_name),"sense"),
    paste(c("All",rep_table$replicate_name),"antisense"),
  )
)

meta_base <- ggplot(mx_data,mapping=aes(x=Position,y=log2FoldChange)) +
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
  ) + 
  scale_colour_manual(
    values=exp_col
  ) +
  scale_fill_manual(
    values=exp_col
  )


meta_plot <- meta_base +
  geom_smooth(
    linewidth=0.8,
    method="loess",
    n=300,
    span=0.05,
    alpha=0.2,
    aes(group=replicate_sense,colour=replicate,linetype=replicate,fill=replicate),
  ) +
  scale_colour_manual(
    breaks=rep_brks,
    values=rep_colours,
  ) +
  scale_fill_manual(
    breaks=rep_brks,
    values=rep_colours,
  ) +
  xlab(
    paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE), "Position (bps)")
  ) +
  scale_x_continuous(
    limits=xlim,
    breaks=meta_xbrks,
  ) +  
  theme(
    legend.position = "right"
  )

if (rep_pair & length(unique(mx_data$rep)) > 1) {

rep_lines <- rep(1,nrow(rep_table))
names(rep_lines) <- rep_table$replicate_name
meta_plot <- meta_plot + 
  facet_wrap(vars(rep),ncol=1,strip.position="top") +
  scale_linetype_manual(
    breaks = rep_brks,
    values = rep_lines 
  )

meta_plot_h <- length(unique(mx_data$rep))*min(c(2.3,7/length(unique(mx_data$rep))))

} else {

meta_plot <- meta_plot +
  scale_linetype_manual(
    breaks = rep_brks,
    values = rep_lines
  )
meta_plot_h <- 3.5
  
}

meta_caption <- paste(
  "Normalised change in log2 ", difference, " coverage over ",meta_gene_n," ", feature_i, " in ", experiment, "." 
, sep="" )
meta_plot_title <- paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE),"change in log2", difference,  "coverage.")
meta_plot_caption <- meta_caption 



mx_data <- mx_df[mx_df$i_group == i,]

meta_gene_n <- sum(rownames(sig_bg)[ (sig_bg$sig2bg >= sig & sig_bg$bg2sig >= bg) ] %in% expr_i$featureID) 


sample_colours <- as.character(unique(mx_data$colour))
sample_names <- as.character(unique(mx_data$Sample))
sample_colours
sample_names
names(sample_colours) <- sample_names
sample_colours

max_val <- ifelse(max(mx_data$value) >= 0, signif(1.1*max(mx_data$value),2), 0)
min_val <- ifelse(min(mx_data$value) >= 0, 0 , signif(1.1*min(mx_data$value),2))
xlim <- c(min(mx_data$Position),max(mx_data$Position))

mx_data$Sample <- gsub("_"," ",mx_data$Sample)

head(mx_data,10)

heat_colours <- c("dodgerblue","blue","black","red","orange")

mx_data$Condition <- factor(mx_data$Condition, levels=c(control, treat))
mx_data$cond_group <- factor(mx_data$cond_group, levels=c(paste(control,"sense"),paste(treat,"sense"),paste(control,"antisense"),paste(treat,"antisense")))

meta <- ggplot(mx_data,mapping=aes(x=Position,y=value,group=cond_group,colour=Condition)) +
  geom_smooth(
    size=0.8,
    method="loess",
    n=300,
    span=0.05,
    alpha=0.2,
    aes(fill=Condition),
  ) +
  scale_colour_manual(
    breaks=c(control,treat),
    values=condition_col,
  ) +
  scale_fill_manual(
    breaks=c(control,treat),
    values=condition_col,
  ) +
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
    panel.border=element_rect(fill=NA,colour="black",size=0.7), 
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",size=0.1),
    legend.background=element_rect(fill="White"), 
    legend.key=element_rect(colour="white",fill="White"), 
    legend.position="none",
    axis.text=element_text(colour="black"),
    axis.line=element_line(colour="black",size=0.1),
    axis.line.x.top=element_line(colour="black",size=0.1),
    axis.line.y.right=element_line(colour="black",size=0.1)
  )



if (section=="body") {

meta <- meta + 
  geom_vline(
    xintercept=main_bin,
    linetype=5,
    colour="black",
    alpha=0.2
  )
}

meta_plot <- meta +
  xlab(
    paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE), "Position (bps)")
  ) +
  scale_x_continuous(
    limits=xlim,
    breaks=meta_xbrks,
  ) +  
  ylab(paste("Normalised",meta_y, "Coverage Depth")) + 
  facet_wrap(vars(rep),ncol=1,strip.position="top") + 
  theme(
    legend.position = "right"
  )
    

meta_caption <- paste(
  "Normalised coverage over ",meta_gene_n," ", feature_i, " in ", experiment, "." 
, sep="" )

meta_plot_title <- paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE),"Coverage Profiles.")
meta_plot_caption <- meta_caption 



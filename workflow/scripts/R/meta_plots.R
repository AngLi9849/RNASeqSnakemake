mx_data <- mx_df[mx_df$i_group == i,]


sample_colours <- as.character(unique(mx_data$colour))
sample_names <- as.character(unique(mx_data$Sample))
sample_colours
sample_names
names(sample_colours) <- sample_names
sample_colours

max_val <- ifelse(max(mx_data$value) >= 0, signif(1.1*max(mx_data$value),2), 0)
min_val <- ifelse(min(mx_data$value) >= 0, 0 , signif(1.1*min(mx_data$value),2))
xlim <- c(min(mx_data$Position),max(mx_data$Position))

if (section=="body") {
xbrks <- c(0,main_bin)
names(xbrks) <- c("Start","End")

xbrk_short <- c(0,main_bin)
names(xbrk_short) <- c("Start","End")

} else {
bef_brk_len <- signif(len_bef_n*1.5,1)/2
bef_brk <- paste("-", as.character(bef_brk_len), sep="")
bef_brk_pos <- floor(signif(bef_bin*1.5,1)/2)


aft_brk_len <- signif(len_bef_n*1.5,1)/2
bef_brk <- paste("+", as.character(aft_brk_len), sep="")
aft_brk_pos <- floor(signif((main_bin-bef_bin)*1.5,1)/2) 


xbrk_short <- c(0)
names(xbrk_short) <- c(toTitleCase(paste(base,section)))

xbrks <- c(ifelse(len_bef_n>0,0-bef_brk_pos,NULL),0,ifelse(len_aft_n>0,aft_brk_pos,NULL))
names(xbrks) <- c(ifelse(len_bef_n>0,bef_brk,NULL),toTitleCase(paste(base,section)),ifelse(len_aft_n>0,aft_brk,NULL))
}
mx_data$Sample <- gsub("_"," ",mx_data$Sample)

head(mx_data,10)


heat_colours <- c("dodgerblue","blue","black","red","orange")

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
    breaks=xbrk_short,
  ) +
  xlab(feature) +
  ylab("Normalised Coverage") +
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
    gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE)
  ) +
  scale_x_continuous(
    limits=xlim,
    breaks=xbrks,
  ) +  
  ylab("Normalised Coverage Depth") + 
  facet_wrap(vars(rep),ncol=1,strip.position="top") + 
  theme(
    legend.position = "right"
  )
    

meta_caption <- paste(
  "Normalised coverage over ", feature_i, " in ", experiment, "." 
, sep="" )

meta_plot_title <- paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE),"Coverage Profiles.")
meta_plot_caption <- meta_caption 



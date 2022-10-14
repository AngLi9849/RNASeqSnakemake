sum_group_n <- length(unique(sum_bar_data$group))

sum_pie_data$Label <- ifelse(sum_pie_data$Numbers > 0,as.character(sum_pie_data$Numbers),NA)

i_group_levels <- unlist(lapply(i_group,function(i) {gsub("exonic ","exonic\n",(gsub("_"," ",toTitleCase(ifelse(i=="","All",i)))))} ) )

sum_pie_data$group <- factor(sum_pie_data$group, levels=i_group_levels)

sum_pie <- ggplot(data = sum_pie_data, aes(x=group, y=Numbers, fill=fct_inorder(Category))) +
  geom_bar(position="fill",stat="identity", colour="black") +
  scale_fill_manual(paste(feature,"\n",toTitleCase(difference),sep=""),values=sum_pie_data$Colours) +
  scale_x_discrete(breaks=sum_pie_data$group) +
  scale_y_continuous(labels=scales::percent) +
  xlab("") +
  ylab("Percentage") +
  geom_text_repel(
    data = sum_pie_data,
    mapping = aes(y=pos, label=fct_inorder(Label)),
    size=2.5,
#    color="black",
    color="white",
    fontface=1,
    bg.color="black",
    bg.r=0.15,
    force=0,
    force_pull=Inf
  ) +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",size=0.1),
    legend.background=element_rect(fill="White"),
    legend.key=element_rect(colour="white",fill="White"),
    axis.text = element_text(colour="black"),
    axis.text.x = if (sum_group_n > 3) (element_text(angle = 45, vjust = 1, hjust=1)) else (element_text()),
    axis.line=element_line(colour="black",size=0.5),
    axis.title.x=element_blank(),
    axis.title.y = element_text(size=9)
  )

bar_cond <- sum_bar_data$condition
names(bar_cond) <- sum_bar_data$condition
bar_max <- (max(sum_bar_data$mean) + max(sum_bar_data$SD))
bar_min <- (min(sum_bar_data$mean) - max(sum_bar_data$SD))
bar_brks <- signif(c(0,1,bar_max),2)

sum_bar_data$Colours <- condition_col[match(bar_cond,names(condition_col))]
sum_bar_data$group <- factor(sum_bar_data$group, levels=i_group_levels)
sum_bar_data$condition <- factor(sum_bar_data$condition,levels=c(control,treat))

sum_bar <- ggplot(data = sum_bar_data, aes(fill=condition, y=mean, x = group)) +
  geom_col(width=0.6,position=position_dodge(0.7),colour = "black") +
  geom_errorbar(aes(ymin=SD_min,ymax=mean+SD),width=0.3,position=position_dodge(0.7)) +
  scale_fill_manual(toTitleCase(gsub("_"," ",protocol)),values=condition_col) + 
  scale_x_discrete(breaks=sum_bar_data$group) +
  scale_y_continuous(limits=c(0,bar_max),breaks=c(0,bar_brks)) +
  geom_hline(
    yintercept=1,
    alpha=0.3,
    linetype=5
  ) +
  ylab(paste("Total",toTitleCase(difference))) +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",size=0.1),
    legend.background=element_rect(fill="White"),
    legend.key=element_rect(colour="white",fill="White"),
    axis.text = element_text(colour="black"),
    axis.text.x = if (sum_group_n > 3) (element_text(angle = 45, vjust = 1, hjust=1)) else (element_text()),
    axis.line=element_line(colour="black",size=0.5),
    axis.title.x=element_blank(),
    axis.title.y = element_text(size=9)
  )

sum_violin_data$group <- factor(sum_violin_data$group, levels=i_group_levels)
sum_violin_data$condition <- factor(sum_violin_data$condition, levels=c(control,treat))


if (difference=="expression levels") {
violin_ymax <- 10^(log10(abs(max(sum_violin_data$value)))^1.2)
violin_p_y <- (abs(max(sum_violin_data$value)))^1.1
} else {
violin_ymax <- max(sum_violin_data$value) + 0.2*abs(max(sum_violin_data$value)-min(sum_violin_data$value))
violin_p_y <- max(sum_violin_data$value) + 0.1*abs(max(sum_violin_data$value)-min(sum_violin_data$value))
}

test_p_i <- compare_means(formula=value ~ condition, data = sum_violin_data, group.by="group", comparisons = compare, method="t.test", paired = TRUE)
test_p_i$y.position <- c(violin_p_y)

sum_violin <- ggplot(data = sum_violin_data, aes(x=group, group=interaction(group,condition), y=value)) +
  geom_violin(trim=TRUE,aes(fill=condition),scale="width",width=0.7) +
  stat_pvalue_manual(data = test_p_i, x="group", label = "p.signif",inherit.aes=F) +
  scale_fill_manual("Conditions",values=condition_col ) +
  geom_boxplot(width=0.1,position=position_dodge(0.7)) +
  scale_x_discrete(breaks=sum_violin_data$group) +
  geom_hline(
    yintercept=1,
    alpha=0.3,
    linetype=5
  ) +
  ylab(paste(toTitleCase(difference))) +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",size=0.1),
    legend.background=element_rect(fill="White"),
    legend.key=element_rect(colour="white",fill="White"),
    axis.text = element_text(colour="black"),
    axis.text.x = if (sum_group_n > 3) (element_text(angle = 45, vjust = 1, hjust=1)) else (element_text()),
    axis.line=element_line(colour="black",size=0.5),
    axis.title.x=element_blank(),
    axis.title.y = element_text(size=9)
  )

if (difference=="expression levels") {
violin_ymin <- min(sum_violin_data$value[!is.infinite(log10(sum_violin_data$value))])
sum_violin <- sum_violin + scale_y_log10(limits = c(violin_ymin, violin_ymax),breaks=c("1x"=1))
} else {
sum_violin <- sum_violin + scale_y_continuous(limits = c(0, violin_ymax),breaks=c("1x"=1))
}
 
overview <- ggarrange(plotlist=list(sum_pie,sum_bar,sum_violin),ncol=1,nrow=3,labels="AUTO")

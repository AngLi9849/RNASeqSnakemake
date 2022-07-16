names(violin_data) <- c("value","condition")
violin_ymax <- max(violin_data$value) + 0.2*abs(max(violin_data$value)-min(violin_data$value))
violin_p_y <- max(violin_data$value) + 0.1*abs(max(violin_data$value)-min(violin_data$value))

violin <- ggplot(data = violin_data, aes(x=condition, y=value)) +
  geom_violin(trim=FALSE,aes(fill=condition)) +
  geom_boxplot(width=0.1) +
  stat_compare_means(comparisons = contrast, label = "p.signif", label.y=violin_p_y, method="wilcox.test", paired = pair) +
  scale_y_continuous(limits = c(0, violin_ymax)) +
  scale_fill_manual("Conditions",values=condition_col[names(condition_col) %in% c(control,exp)], labels = bar_cond) +
  ylab(paste(feature,change)) +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",size=0.1),
    legend.background=element_rect(fill="White"),
    legend.key=element_rect(colour="white",fill="White"),
    legend.position="none",
    axis.text = element_text(colour="black"),
    axis.line=element_line(colour="black",size=0.5),
    axis.title.x=element_blank(),
    axis.title.y = element_text(size=9)
  )


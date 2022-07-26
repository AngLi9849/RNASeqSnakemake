violin_data <- data.frame(
  unlist(
    lapply(colnames(mean_level_i), function (x) {
      mean_level_i[paste(x)]
    })),
  unlist(
    lapply(colnames(mean_level_i), function (x) {
      replicate(nrow(mean_level_i),paste(x))
    }))
)

names(violin_data) <- c("value","condition")

if (difference=="expression levels") {
violin_ymax <- max(violin_data$value) + (10^(0.2)-1)*(abs(max(violin_data$value)-min(violin_data$value)))
violin_p_y <- max(violin_data$value) + (10^(0.1)-1)*abs(max(violin_data$value)-min(violin_data$value))
} else {
violin_ymax <- max(violin_data$value) + 0.2*abs(max(violin_data$value)-min(violin_data$value))
violin_p_y <- max(violin_data$value) + 0.1*abs(max(violin_data$value)-min(violin_data$value))
}

violin <- ggplot(data = violin_data, aes(x=condition, y=value)) +
  geom_violin(trim=FALSE,aes(fill=condition)) +
  geom_boxplot(width=0.1) +
  stat_compare_means(comparisons = compare, label = "p.signif", label.y=violin_p_y, method="t.test", paired = TRUE) +
  scale_fill_manual("Conditions",values=condition_col[names(condition_col) %in% c(control,treat)], labels = c(control,treat) ) +
  ylab(paste("Normalised",toTitleCase(difference))) +
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

if (difference=="expression levels") {
violin <- violin + scale_y_log10(limits = c(1, violin_ymax))
} else {
violin <- violin + scale_y_continuous(limits = c(0, violin_ymax))
}

test_p_i <- compare_means(formula=value ~ condition, data = violin_data, comparisons = compare, method="t.test", paired = TRUE)

violin_caption <- paste(
  "Distribution of ", difference, " of ", nrow(cts_genes_i), " ",  feature_i, " in violin and bar plot. Significance of changes by ", as.character(test_p_i$method[1]), " method is indicated. (p=", test_p_i$p.format[1] , ", ", test_p_i$p.signif[1],")", sep=""
  )

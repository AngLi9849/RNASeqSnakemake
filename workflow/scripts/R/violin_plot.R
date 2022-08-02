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

if (difference=="expression") {
violin_ymax <- 10^(log10(abs(max(violin_data$value)))^1.2)
violin_p_y <- log10(abs(max(violin_data$value)))^1.1
} else {
violin_ymax <- max(violin_data$value) + 0.2*abs(max(violin_data$value)-min(violin_data$value))
violin_p_y <- max(violin_data$value) + 0.1*abs(max(violin_data$value)-min(violin_data$value))
}

test_p_i <- compare_means(formula=value ~ condition, data = violin_data, comparisons = compare, method="t.test", paired = TRUE)
test_p_i$y.position <- c(violin_p_y)

violin <- ggplot(data = violin_data, aes(x=condition, y=value)) +
  geom_violin(trim=FALSE,aes(fill=condition)) +
  geom_boxplot(width=0.1) +
  stat_pvalue_manual(data = test_p_i, label = "p.signif") +
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

if (difference=="expression") {
violin_ymin <- min(violin_data$value[!is.infinite(log10(violin_data$value))])
violin <- violin + scale_y_log10(limits = c(1, violin_ymax))
} else {
violin <- violin + scale_y_continuous(limits = c(0, violin_ymax))
}

test_p_i <- compare_means(formula=value ~ condition, data = violin_data, comparisons = compare, method="t.test", paired = TRUE)

violin_caption <- paste(
  "Distribution of normalised ", difference, " of ", total_i, " ",  feature_i, " presented as violin and bar plot. Indicated significance of comparison is computed with ", as.character(test_p_i$method[1]), " method. (p=", test_p_i$p.format[1] , ", ", test_p_i$p.signif[1],")", sep=""
  )

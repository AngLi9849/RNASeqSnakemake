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

violin_data$condition <- factor(violin_data$condition, levels=c(control,treat))

if (difference=="expression levels") {
violin_ymax <- 10^(log10(abs(max(violin_data$value)))^1.2)
violin_p_y <- log10(abs(max(violin_data$value)))^1.1
} else {
violin_ymax <- max(violin_data$value) + 0.2*abs(max(violin_data$value)-min(violin_data$value))
violin_p_y <- max(violin_data$value) + 0.1*abs(max(violin_data$value)-min(violin_data$value))
}

violin_data_i <- violin_data
violin_data_i$group <- group_label

violin_control_mean <- median(violin_data_i$value[violin_data_i$condition==control])
violin_data_i$value <- violin_data_i$value/violin_control_mean 
violin_data_i$Colours <- condition_col[match(violin_data_i$condition,names(condition_col))]

if (nrow(mean_level_i) > 5) {
  if (exists("sum_violin_data")) {

  sum_violin_data <- rbind(sum_violin_data,violin_data_i)

  } else {

  sum_violin_data <- violin_data_i

  }
}

violin <- ggplot(data = violin_data, aes(x=condition, y=value)) +
  geom_violin(trim=FALSE,aes(fill=condition)) +
  geom_boxplot(width=0.1) +
  scale_fill_manual("Conditions",values=condition_col[names(condition_col) %in% c(control,treat)], labels = c(control,treat) ) +
  ylab(paste("Mean",toTitleCase(difference_unit))) +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",linewidth=0.1),
    legend.background=element_rect(fill="White"),
    legend.key=element_rect(colour="white",fill="White"),
    legend.position="none",
    axis.text = element_text(colour="black"),
    axis.line=element_line(colour="black",linewidth=0.5),
    axis.title.x=element_blank(),
    axis.title.y = element_text(size=9)
  )

if (difference=="expression levels") {
violin_ymin <- min(violin_data$value[!is.infinite(log10(violin_data$value))])
violin <- violin + scale_y_log10(limits = c(violin_ymin, violin_ymax))
} else {
violin <- violin + scale_y_continuous(limits = c(0, violin_ymax))
}

violin_caption <- paste(
  "Distribution of normalised ", difference, " of ", nrow(mean_level_i), " expressed ",  feature_i, " presented as violin and box plot.", 
sep="") 

if (nrow(mean_level_i) > 20) {
  test_p_i <- compare_means(formula=value ~ condition, data = violin_data, comparisons = compare, method="t.test", paired = TRUE)
  test_p_i$y.position <- c(violin_p_y)


violin <- violin + stat_pvalue_manual(data = test_p_i, label = "p.signif")

violin_caption <- paste( 
  violin_caption, 
  "Indicated significance of comparison is computed with ", as.character(test_p_i$method[1]), " method. (p=", test_p_i$p.format[1] , ", ", test_p_i$p.signif[1],")", sep=""
)

}


dot_data <- data.frame(
  unlist(
    lapply(colnames(mean_level_i), function (x) { 
      mean_level_i[paste(x)] 
    })),
  unlist(
    lapply(colnames(mean_level_i), function (x) { 
      replicate(nrow(mean_level_i),paste(x)) 
    })),
  unlist(
    lapply(colnames(mean_level_i), function (x) { 
      rownames(mean_level_i) 
    }))
)

names(dot_data) <- c("value","condition","featureID")
dot_data$rpkm <- expr$RPKM[match(dot_data$featureID,expr$featureID)]
dot_data$top10 <- ifelse(dot_data$rpkm>=quantile(dot_data$rpkm,0.9),dot_data$value,NA)
dot_data$top50 <- ifelse(dot_data$rpkm>=quantile(dot_data$rpkm,0.5),dot_data$value,NA)

if (difference=="expression levels") {
dot_ymax <- 10^(log10(abs(max(dot_data$value)))^1.2)
dot_p_y <- log10(abs(max(dot_data$value)))^1.1
} else {
dot_ymax <- max(dot_data$value) + 0.2*abs(max(dot_data$value)-min(dot_data$value))
dot_p_y <- max(dot_data$value) + 0.1*abs(max(dot_data$value)-min(dot_data$value))
}


dot_p_y <- max(dot_data$value) + 0.1*abs(max(dot_data$value)-min(dot_data$value))
dot_data$condition <- factor(dot_data$condition, levels=c(control,treat))

dot_bin <- abs(max(dot_data$value) - min(dot_data$value))/(max(c(nrow(dot_data),4000))/50)

dot <- ggplot(data = dot_data, aes(x=factor(condition,levels=c(control,treat)),y=value,fill=condition,colour=condition)) +
  geom_dotplot(
    aes(
      x=dot_data$condition,
      y=dot_data$value
    ),
    binwidth=dot_bin,
    binaxis="y",
    stackdir="center",
    method="histodot",
    alpha=0.05
  ) +
  geom_dotplot(
    aes(
      x=dot_data$condition,
      y=dot_data$top50
    ),
    binwidth=dot_bin,
    binaxis="y",
    stackdir="center",
    method="histodot",
    alpha=0.1
  ) +
  geom_dotplot(
    aes(
      x=dot_data$condition,
      y=dot_data$top10
    ),
    binwidth=dot_bin,
    binaxis="y",
    stackdir="center",
    method="histodot",
    alpha=1
  ) +
  scale_fill_manual("Conditions",values=condition_col[names(condition_col) %in% c(control,treat)], labels = c(control,treat)) +
  scale_colour_manual("Conditions",values=condition_col[names(condition_col) %in% c(control,treat)], labels = c(control,treat)) +
  ylab(paste("Mean",toTitleCase(difference_unit))) +
  stat_compare_means(comparisons = contrast, label = "p.signif", label.y=dot_p_y, method="t.test") +
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
dot_ymin <- min(dot_data$value[!is.infinite(log10(dot_data$value))])
dot <- dot + scale_y_log10(limits = c(dot_ymin, dot_ymax))
} else {
dot <- dot + scale_y_continuous(limits = c(0, dot_ymax))
}

dot_caption <- paste(
  "Distribution of normalised ", difference, " of ", nrow(mean_level_i), " expressed ",  feature_i, " in dot plot. Top 50 and 10 percent highliest expressed ", feature_i, " are represented by increasingly opaque dots.", sep=""
  )
 
#ggsave("dot.png",dot,width=w/2,height=h/4,units="cm",dpi=300)


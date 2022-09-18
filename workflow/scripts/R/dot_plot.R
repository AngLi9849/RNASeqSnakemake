
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
dot_data$rpkm <- expr$rpkm[match(dot_data$featureID,expr$featureID)]
dot_data$top10 <- ifelse(dot_data$rpkm>=quantile(dot_data$rpkm,0.9),dot_data$value,NA)
dot_data$top50 <- ifelse(dot_data$rpkm>=quantile(dot_data$rpkm,0.5),dot_data$value,NA)


dot_ymax <- max(dot_data$value) + 0.2*abs(max(dot_data$value)-min(dot_data$value))
dot_ymin <- min(dot_data$value)
dot_p_y <- max(dot_data$value) + 0.1*abs(max(dot_data$value)-min(dot_data$value))
dot_data$condition <- factor(dot_data$condition, levels=c(control,treat))

dot_bin <- abs(max(dot_data$value) - min(dot_data$value))/(nrow(dot_data)/50)
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
  ylab("Splicing Ratio") +
  stat_compare_means(comparisons = contrast, label = "p.signif", label.y=dot_p_y, method="t.test") +
  scale_y_continuous(limits = c(dot_ymin, dot_ymax)) +
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

dot_caption <- paste(
  "Distribution of normalised ", difference, " of ", nrow(mean_level_i), " expressed ",  feature_i, " in dot plot. Top 50 and 10 percent highliest expressed ", feature_i, " are represented by increasingly opaque dots.", sep=""
  )
 
#ggsave("dot.png",dot,width=w/2,height=h/4,units="cm",dpi=300)


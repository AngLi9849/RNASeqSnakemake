bias_bin <- 100
bias_lim <- c(1,bias_bin)
sample_colours <- sample_table$colour
names(sample_colours) <- gsub("_"," ",sample_table$sample_name)
sample_brks <- gsub("_"," ",sample_table$sample_name)
biases_n <- 0
for (b in biases$bias) {

biases_n <- biases_n + 1

expr_full$ntile <- ntile(x=unlist(expr_full[b]),n=bias_bin)
count_data <- cts_i
count_data['rank'] <- expr_full$ntile[match(rownames(count_data),expr_full$featureID)]
count_data <- aggregate(count_data[colnames(count_data)!="rank"],by=list(rank=count_data$rank),FUN=sum)

count_sum <- count_data[colnames(count_data)!="rank"]

count_bias_data <- data.frame(
  unlist(
    lapply(colnames(count_sum), function (x) {
      count_data[paste(x)]
    })),
  unlist(
    lapply(colnames(count_sum), function (x) {
      replicate(nrow(count_sum),gsub("_"," ",x))
    })),
  rep(count_data$rank,ncol(count_sum))
)


names(count_bias_data) <- c("Counts","Sample","Rank")

count_bias_data$replicate <- sample_table$replicate[match(count_bias_data$Sample,gsub("_"," ",sample_table$sample_name))]

b_unit <- biases$unit[biases$bias==b]
bias_lab <- paste(b, ifelse(b_unit=="", "",paste("(", b_unit ,")",sep="")),sep="")

bias_brks <- c(
  1,
  bias_bin/4,
  bias_bin/2,
  bias_bin*0.75,
  bias_bin
)

names(bias_brks) <- paste(
  signif(
    c(
      min(expr_full[b],na.rm=T),
      quantile(expr_full[b],0.25,na.rm=T),
      quantile(expr_full[b],0.5,na.rm=T),
      quantile(expr_full[b],0.75,na.rm=T),
      max(expr_full[b],na.rm=T)
    ) * ifelse(as.logical(biases$percent[biases$bias==b]), 100, 1),
    3
  ), 
  ifelse(as.logical(biases$percent[biases$bias==b]), "%", ""), 
  sep=""
)

sample_lines <- sample_table$replicate
names(sample_lines) <- gsub("_"," ",sample_table$sample_name)

sample

count_bias <- ggplot(data=count_bias_data, aes(x=Rank,y=Counts)) +
  geom_smooth(
    linewidth=0.8,
    method="loess",
    n=300,
    span=0.5,
    alpha=0.2,
    aes(colour=Sample,linetype=Sample,fill=Sample)
  ) +
  scale_colour_manual(
    breaks=sample_brks,
    values=sample_colours,
  ) +
  scale_fill_manual(
    breaks=sample_brks,
    values=sample_colours,
  ) +
  scale_x_continuous(
    limits=bias_lim,
    breaks=bias_brks,
  ) +
  scale_linetype_manual(
    breaks = sample_brks,
    values = sample_lines,
  ) +
  labs(
    subtitle=paste(difference)
  ) +
  xlab(
    paste(feature_i,bias_lab)
  ) +
  ylab("Read Counts") +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    panel.border=element_rect(fill=NA,colour="black",linewidth=0.7),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",linewidth=0.1),
    legend.background=element_rect(fill="White"),
    legend.key=element_rect(colour="white",fill="White"),
    legend.position=ifelse(biases_n==1,"bottom","none"),
    legend.direction="vertical",
    legend.just="left",
    axis.text=element_text(colour="black"),
    axis.text.x = if (length(bias_brks)>6) (element_text(angle = 45, vjust = 1, hjust=1,colour="black")) else (element_text(colour="black")),
    axis.line=element_line(colour="black",linewidth=0.1),
    axis.line.x.top=element_line(colour="black",linewidth=0.1),
    axis.line.y.right=element_line(colour="black",linewidth=0.1)
  )

if (biases_n == 1) {
  count_bias_legend <- as_ggplot(get_legend(count_bias))
  count_bias <- count_bias + 
    theme(
      legend.position="none"
    )
}

assign(paste(b,"_count_bias",sep=""),count_bias)

}

count_bias_ls <- lapply(c(paste(biases$bias,"_count_bias",sep=""),"count_bias_legend"),function(x){get(x)})

count_bias <- ggarrange(plotlist=count_bias_ls,ncol=1,nrow=nrow(biases)+1)

count_bias_title <- "Count Bias"
count_bias_caption <- analysis_title
count_bias_h <- 7.5


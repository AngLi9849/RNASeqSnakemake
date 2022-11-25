bias_bin <- 100
bias_lim <- c(1,bias_bin)
sample_colours <- sample_table$colour
names(sample_colours) <- sample_table$sample_name
sample_brks <- gsub("_"," ",sample_table$sample_name)

for (bias in bias$bias) {

expr_full$ntile <- ntile(x=expr_full$bias,n=bias_bin)
count_data <- cts_i
count_data['rank'] <- expr_full$ntile[match(rownames(count_data),expr_full$featureID)]
count_data <- aggregate(count_data[colnames(count_data)!=rank],by=list(rank=count_data$rank),FUN=sum)

count_sum <- count_data[colnames(count_data)!="rank"]

count_bias_data <- data.frame(
  unlist(
    lapply(colnames(count_sum), function (x) {
      count_data[paste(x)]
    })),
  unlist(
    lapply(colnames(count_sum), function (x) {
      replicate(nrow(count_sum),paste(x))
    })),
  rep(count_data$rank,ncol(count_data))
)


names(count_bias_data) <- c("Counts","Sample","Rank")

count_bias_data$replicate <- sample_table$replicate[match(count_bias_data$Sample,sample_table$sample_name)]

bias_lab <- paste(bias, " (", bias$unit[bias$bias==bias],")",sep="")

bias_brks <- c(
  1,
  bias_bin/4,
  bias_bin/2,
  bias_bin*0.75,
  bias_bin
)

names(bias_brks) <- c(
  min(expr_full[bias],na.rm=T),
  quantile(expr_full[bias],0.25,na.rm=T),
  quantile(expr_full[bias],0.5,na.rm=T),
  quantile(expr_full[bias],0.75,na.rm=T),
  max(expr_full[bias],na.rm=T),
}

count_bias <- ggplot(data=count_bias_data, aes(x=Rank,y=Counts,colour=Sample,linetype=replicate))
  geom_smooth(
    size=0.8,
    method="loess",
    n=300,
    span=0.05,
    alpha=0.2,
    aes(fill=Sample),
  ) +
  scale_colour_manual(
    breaks=sample_brks,
    values=sample_colours,
  ) +
  scale_x_continuous(
    limits=bias_lim,
    breaks=bias_brks,
  ) +
  labs(
    subtitle=paste()
  ) +
  xlab("Read Counts") +
  ylab(bias_lab) +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    panel.border=element_rect(fill=NA,colour="black",size=0.7),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",size=0.1),
    legend.background=element_rect(fill="White"),
    legend.key=element_rect(colour="white",fill="White"),
    legend.position="none",
    axis.text=element_text(colour="black"),
    axis.text.x = if (length(heat_xbrks)>2) (element_text(angle = 45, vjust = 1, hjust=1,colour="black")) else (element_text(colour="black")),
    axis.line=element_line(colour="black",size=0.1),
    axis.line.x.top=element_line(colour="black",size=0.1),
    axis.line.y.right=element_line(colour="black",size=0.1)
  )

assign(paste(bias,"_count_bias",sep=""),count_bias)
assign(paste(bias,"_count_bias_titile",sep=""),paste(bias,"Count Bias"))
bias_caption <- ""
assign(paste(bias,"_count_bias_caption",sep=""),bias_caption)

}

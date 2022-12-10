bias_bin <- 100
bias_lim <- c(1,bias_bin)
sample_colours <- sample_table$colour
names(sample_colours) <- gsub("_"," ",sample_table$sample_name)
sample_brks <- gsub("_"," ",sample_table$sample_name)

for (b in biases$bias) {

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

bias_lab <- paste(bias, " (", biases$unit[biases$bias==b],")",sep="")

bias_brks <- c(
  1,
  bias_bin/4,
  bias_bin/2,
  bias_bin*0.75,
  bias_bin
)

names(bias_brks) <- c(
  min(expr_full[b],na.rm=T),
  quantile(expr_full[b],0.25,na.rm=T),
  quantile(expr_full[b],0.5,na.rm=T),
  quantile(expr_full[b],0.75,na.rm=T),
  max(expr_full[b],na.rm=T)
)

sample_lines <- sample_table$replicate
names(sample_lines) <- gsub("_"," ",sample_table$sample_name)

sample

count_bias <- ggplot(data=count_bias_data, aes(x=Rank,y=Counts)) +
  geom_smooth(
    size=0.8,
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
  xlab(bias_lab) +
  ylab("Read Counts") +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    panel.border=element_rect(fill=NA,colour="black",size=0.7),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",size=0.1),
    legend.background=element_rect(fill="White"),
    legend.key=element_rect(colour="white",fill="White"),
    legend.position="bottom",
    legend.direction="vertical",
    legend.just="left",
    axis.text=element_text(colour="black"),
    axis.text.x = if (length(heat_xbrks)>2) (element_text(angle = 45, vjust = 1, hjust=1,colour="black")) else (element_text(colour="black")),
    axis.line=element_line(colour="black",size=0.1),
    axis.line.x.top=element_line(colour="black",size=0.1),
    axis.line.y.right=element_line(colour="black",size=0.1)
  )

assign(paste(b,"_count_bias",sep=""),count_bias)
assign(paste(b,"_count_bias_title",sep=""),paste(b,"Count Bias"))
bias_caption <- analysis_title
assign(paste(b,"_count_bias_caption",sep=""),bias_caption)
assign(paste(b,"_count_bias_h",sep=""),7.5)

}

sig_up <- sum(expr_i$log2FoldChange[expr_i$padj < sig_p] > 0,na.rm=T)
sig_down <- sum(expr_i$log2FoldChange[expr_i$padj < sig_p] < 0,na.rm=T)
insig <- sum(expr_i$padj[expr_i$padj < undetect_p] >= sig_p,na.rm=T)
undetect <- sum(expr_i$padj >= undetect_p,na.rm=T) 
insuf_i <- max(total_i - undetect - sig_up - sig_down - insig,0)
total_i <- total_i - insuf_i


sig_up
sig_down
insig
undetect

sig_up_pc <- paste(signif(sig_up/total_i*100,2), "%",sep="")
sig_down_pc <- paste(signif(sig_down/total_i*100,2), "%",sep="")
insig_pc <- paste(signif(insig/total_i*100,2), "%",sep="")
undetect_pc <- paste(signif(undetect/total_i*100,2), "%",sep="")
insuf_pc <- paste(signif(insuf_i/total_i*100,2), "%",sep="")

pie_data <- data.frame(
  data1=c("Insufficient Reads", "Unchanged","Insignificant","Increase","Decrease"),
  data2=c(insuf_i,undetect,insig,sig_up,sig_down),
  data3=c(insuf_pc,undetect_pc,insig_pc,sig_up_pc,sig_down_pc),
  data3=c(void_col,background_col,insig_col,up_col,down_col)
)



names(pie_data) <- c("Category","Numbers","Percent","Colours")
pie_label <- paste(total_i,"Evidenced")
pie_data$Percent <- paste(signif(pie_data$Numbers/total_i*100,2), "%")
pie_data$Label <- paste(pie_data$Numbers,pie_data$Category)

pie_data <- pie_data[pie_data$Category != "Insufficient Reads",]



pie_data <- pie_data %>%
  mutate(csum = rev(cumsum(rev(Numbers))),
         pos = Numbers/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Numbers/2, pos))

pie_data_i <- pie_data
pie_data_i$group <- group_label
pie_data_i$pos <- pie_data_i$pos/total_i

if (exists("sum_pie_data")) {

sum_pie_data <- rbind(sum_pie_data,pie_data_i)

} else {

sum_pie_data <- pie_data_i

}


pie <- ggplot(data = pie_data, aes(x="", y=Numbers, fill=fct_inorder(Label))) +
  geom_bar(width=1,stat="identity", colour="black") +
  coord_polar("y",start=0) +
  scale_fill_manual(paste(feature_i),values=pie_data$Colours) +
  xlab("") +
  ylab(pie_label) +
  geom_label_repel(
    data = pie_data,
    mapping = aes(y=pos, label=fct_inorder(Label)),
    size=3,
    color="black",
    fontface=1,
    segment.color="black",
    segment.alpha=1,
    label.size = NA,
    fill = alpha(c("white"),0.8),
    nudge_x=0.7,
    force = 10,
  ) +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",linewidth=0.1),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line=element_blank(),
    axis.title.x = element_text(size=9)
  )

rpkm_lim <- ifelse(round(min_rpkm,3)==0,").",paste(" or ",round(min_rpkm,3), " RPKM). ",sep=""))

pie_caption <- paste(
  "In ", total_i, " ", feature_i, ", ", difference, " of ",
  sig_up, " significantly increased (", sig_up_pc, ", ", up_col, "), ", 
  sig_down, " significantly decreased (", sig_down_pc, ", ", down_col, "), and ",
  insig, " changed insignificantly (", insig_pc, ", p >= ", sig_p, "). ",
  "Changes were not detected in ", undetect, " (", undetect_pc, ", p >= ", undetect_p, "), and ",
  insuf_i, " were insufficiently evidenced (", "counted less than ", round(min_mean,2), " reads per sample", rpkm_lim, 
  sep="")



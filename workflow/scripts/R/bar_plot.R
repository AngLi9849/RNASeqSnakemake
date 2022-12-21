bar_data <- data.frame(c(control,treat),
                      as.numeric(paste(lapply(c(control_cond,treatment),function(x) {mean(sum_i[rownames(sum_i) %in% sample_table$sample_name[sample_table$condition==x],])}))),
                      as.numeric(paste(lapply(c(control_cond,treatment),function(x) {sd(sum_i[rownames(sum_i) %in% sample_table$sample_name[sample_table$condition==x],])})))
)

names(bar_data) <- c("condition","mean","SD")
bar_cond <- bar_data$condition
names(bar_cond) <- bar_data$condition
bar_max <- (max(bar_data$mean) + max(bar_data$SD))
bar_min <- (min(bar_data$mean) - max(bar_data$SD))
bar_brks <- signif(c(bar_max/2,bar_max),2)

bar_data$SD_min <- ifelse(bar_data$SD > bar_data$mean, bar_data$mean, bar_data$mean-bar_data$SD)

bar_data_i <- bar_data
bar_data_i$SD_pc <- (bar_data_i$SD/bar_data_i$mean)
bar_data_i$mean <- bar_data_i$mean/(bar_data_i$mean[bar_data_i$condition==control])
bar_data_i$SD <- bar_data_i$mean*bar_data_i$SD_pc
bar_data_i$SD_min <- ifelse(bar_data_i$SD > bar_data_i$mean, bar_data_i$mean, bar_data_i$mean-bar_data_i$SD)
bar_data_i$group <- group_label

if (exists("sum_bar_data")) {

sum_bar_data <- rbind(sum_bar_data,bar_data_i)

} else {

sum_bar_data <- bar_data_i

}

bar <- ggplot(data = bar_data, aes(x=factor(condition,levels=c(control,treat)), y=mean, fill = condition)) + 
  geom_col(width=0.6,colour = "black") + 
  geom_errorbar(aes(ymin=SD_min,ymax=mean+SD),width=0.2,position=position_dodge()) +
  scale_fill_manual("Conditions",values=condition_col[names(condition_col) %in% c(control,treat)], labels = bar_cond) +
  scale_x_discrete(breaks=bar_data$condition) +
  scale_y_continuous(limits=c(0,bar_max),breaks=c(0,bar_brks)) +
  ylab(paste("Total",toTitleCase(difference_unit))) +
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

bar_caption <- paste(
  "Bar chart showing total ", feature_i, " ", difference, " in ", experiment, ". Error bars represent standard deviations.", sep=""
) 

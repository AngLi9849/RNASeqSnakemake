pie_data <- pie_data %>%
  mutate(csum = rev(cumsum(rev(Numbers))),
         pos = Numbers/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Numbers/2, pos))


pie <- ggplot(data = pie_data, aes(x="", y=Numbers, fill=fct_inorder(Category))) +
  geom_bar(width=1,stat="identity", colour="black") +
  coord_polar("y",start=0) +
  scale_fill_manual(paste(feature,change),values=pie_data$Colours) +
  xlab("") +
  ylab(pie_label) +
  geom_label_repel(
    data = pie_data,
    mapping = aes(y=pos, label=fct_inorder(Category)),
    size=2,
    color="black",
    fontface=1,
    segment.color="black",
    segment.alpha=1,
    label.size = NA,
    fill = alpha(c("white"),0.8),
    nudge_x=0.5,
    force = 10,
  ) +
  theme(
    panel.background=element_rect(fill="White",colour="white"),
    strip.text=element_text(face="bold"),
    strip.background=element_rect(colour="white",fill="white",size=0.1),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line=element_blank(),
    axis.title.y = element_text(size=9)
  )


heat_data <- bind_rows(lapply(lfc_i,function(x) {data.frame(x[,c("heat","title","rank")])}))
names(heat_data) <- c("heat","title","rank")

heat_data$title <- factor(heat_data$title,levels=titles)

heat_scale <- quantile(abs(heat_data$heat),heat_scale_pc,na.rm=T)
heat_data$heat <- heat_data$heat/heat_scale
heat_data$heat <- ifelse(heat_data$heat >= 1, 1, ifelse(heat_data$heat <= -1, -1, heat_data$heat))

heat_lfcbrks_i <- heat_lfcbrks/heat_scale

heat_name <- paste("n=",gene_n,"\nlog2FC",sep="")

heatmap <- ggplot(
  heat_data,
  aes(x=title,y=rank,fill=heat)
  ) +
geom_raster() +
scale_fill_gradientn(
  name=heat_name,
  colours = heat_colours,
  limits=c(-1,1),
  breaks=heat_lfcbrks_i,
  guide = guide_colorbar(
    label = TRUE,
    draw.ulim = TRUE,
    draw.llim = TRUE,
    frame.colour = "black",
    frame.linewidth = 0.3,
    ticks.colour = "black",
    ticks.linewidth = 0.3,
    ticks = TRUE,
    label.position = "right",
    barwidth = 0.5,
    barheight = 6,
    direction = 'vertical')
  ) +
scale_y_continuous(
  limits=c(0.5,max(heat_data$rank)+0.5),
  breaks=c(0),
  expand=c(0,0)
  ) +
scale_x_discrete(
  expand=c(0,0),
  breaks=titles,
  position="top"
  ) +
xlab(paste(group_title,"\n",i_heading," genes")) +
ylab(heat_ylab) +
theme(
  panel.background=element_rect(
    fill="White",
    colour="white"
    ),
  panel.border=element_rect(
    fill=NA,
    colour="black",
    linewidth=0.7
    ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.background=element_rect(
    fill="White"
    ),
  legend.key=element_rect(
    colour="white",
    fill="white"
    ),
  axis.line=element_line(
    colour=NA,
    linewidth=NA
    ),
  axis.title.y = element_text(size=9),
  axis.title.x = element_text(size=9),
  legend.title = element_text(size=9),
  strip.text=element_text(face="bold"),
  strip.background=element_rect(colour="white",fill="white",linewidth=0.1),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, vjust = 0, hjust=0,colour="black")
  )


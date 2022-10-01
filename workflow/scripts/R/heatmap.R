heat_data <- head_df[heat_df$featureID %in% expr_i$featureID,]

heat_ranks <- c("Fold_Change","GC","Length","RPKM")





for (i in heat_ranks) {

heat_data$Rank <- c(expr_i[,paste(i,"_Rank",sep='\t')])[match(heat_data$featureID,expr$featureID),] 



heatmap <- ggplot(
  heat_data,
  aes(x=Position,y=Rank,fill=heat)
  ) +
geom_raster() +
scale_fill_gradientn(
  colours = heat_colours,
  limits=c(-1,1),
  breaks=heat_lfcbrks,
  guide = guide_colorbar(
    label = TRUE,
    draw.ulim = TRUE, 
    draw.llim = TRUE,
    frame.colour = "black",
    frame.linewidth = 1,
    ticks.colour = "black",
    ticks.linewidth = 1,
    ticks = TRUE, 
    label.position = "right",
    barwidth = 0.5,
    barheight = 6,
    direction = 'vertical')
  ) +
scale_y_continuous(
  limits=c(0.5,max(heat_data$Rank)+0.5),
  breaks=NULL,
  expand=c(0,0)
  ) +
scale_x_continuous(
  limits=xlim,
  breaks=xbrks,
  expand=c(0,0)
  ) +
xlab("Relative Position (bps)") +
ylab("Total log2 Fold Change") +
theme(
  panel.background=element_rect(
    fill="White",
    colour="white"
  ),
  panel.border=element_rect(
    fill=NA,
    colour="black",
    size=0.7
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
    size=NA)
  )

}



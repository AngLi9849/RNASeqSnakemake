heat_data <- head_df[heat_df$featureID %in% expr_i$featureID,]

heat_ranks <- c("Fold_Change","GC","Length","RPKM")

heat_ylab <- data.frame(heat_ranks, c(paste(difference, "Fold Change (Low to High)"), "GC Content (%)", paste(title_feature_i, "Mean Expression Levels (RPKM)",sep=" "),paste(ifelse(use_base_length,title_base_i,title_feature_i), "Length (bps)",sep=" ")))

heat_ls <- paste(heat_ranks,"_heat",sep="")

names(heat_ylab) <- c("Ranking","ylab") 

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
  limits=heat_xlim,
  breaks=heat_xbrks,
  expand=c(0,0)
  ) +
xlab(
  paste(gsub("(?<!\\w)(.)","\\U\\1", feature_i, perl = TRUE), "Position (bps)")
  ) +
ylab(heat_ylab$ylab[heat_ylab$Ranking==i]) +
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

assign(paste(i,"_heat",sep=""),heatmap) 

}

heat_ls <- lapply(heat_ls,get)

heat_plot <- ggarrange(plotlist=heat_ls,ncol=2,nrow=lenth(heat_ls)/2,labels="AUTO")

heat_plot_caption <- paste(
  "Heatmaps of fold changes in ", difference, " of ", feature_i, " in ", experiment, ".",
  sep="")



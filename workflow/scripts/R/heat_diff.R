heat_config <- biases[as.logical(biases$heat),]


for (r in heat_config$bias) {
expr_heat <- expr_heat %>% arrange(!!sym(r)) %>% mutate(!!paste(r,"_Rank",sep='') := 1:n()) %>% ungroup
}


heat_pc = c(0,25,50,75,100)
heat_ls <- paste(heat_config$bias,"_heat",sep="")
heat_gene_n <- length(unique(heat_data$featureID))

heat_name <- paste("log2FC")

heat_data$Sense <- factor(heat_data$Sense, levels=c("Sense","Antisense"))

heat_scale <- max(quantile(abs(heat_data$heat[heat_data$heat!=0]),heat_scale_pc,na.rm=T)/1,heat_median_pc/median(abs(heat_data$heat[heat_data$heat!=0 & abs(heat_data$heat)!=1]),na.rm=T))
heat_data$heat <- heat_data$heat*heat_scale
heat_data$heat <- ifelse(heat_data$heat >= 1, 1, ifelse(heat_data$heat <= -1, -1, heat_data$heat))

heat_lfcbrks_i <- heat_lfcbrks*heat_scale
#heat_data$mx_heat <- (((heat_data$heat+1)*heat_data$coverage) +1)/(heat_data$coverage +1 )-1
if(length(heat_lfcbrks_i[abs(heat_lfcbrks_i)<=min(abs(heat_lfcbrks_val[heat_lfcbrks_val!=0]))]) < 3) {
  heat_lfcbrks_temp <- signif((log2((1+0.9/heat_scale)/(1-(0.9/heat_scale)))*2),1)/2
  heat_lfcbrks_temp <- c(0-heat_lfcbrks_temp,0,heat_lfcbrks_temp)
  heat_lfcbrks_i <- unlist(lapply(heat_lfcbrks_temp, function(x) {((2^(x+1))/((2^x)+1))-1}))
  heat_lfcbrks_i <- heat_lfcbrks_i * heat_scale
  names(heat_lfcbrks_i) <- heat_lfcbrks_temp
}

for (r in heat_config$bias) {
heat_unit <- heat_config$unit[heat_config$bias==r]
r_unit <- ifelse(heat_unit=="","",paste("(",heat_unit,")",sep=""))
heat_data$Rank <- as.numeric(unlist(expr_heat[match(heat_data$featureID,expr_heat$featureID),paste(r,"_Rank",sep='')]))
heat_ybrks <- unlist(lapply(heat_pc, function(p) { quantile(as.numeric(unlist(expr_heat[,paste(r, "_Rank",sep='')])),p/100,na.rm=T) } ) )
names(heat_ybrks) <-  unlist(lapply(heat_pc, function(p) {ifelse(heat_unit=="%",100,1)*signif(quantile(as.numeric(unlist(expr_heat[,paste(r)])),p/100,na.rm=T),2) } ) )

heat_gene_lab <- heat_data$Rank
names(heat_gene_lab) <- heat_data$root_name

r_lab <- ifelse(r =="log2FoldChange", paste(toTitleCase(difference), "log2FoldChange"),r)
r_lab <- ifelse(r_lab =="Length", paste(ifelse(use_base_length,title_base,title_feature),"Length"),r_lab)

heat_ylab <- paste(gsub("_"," ",gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",r_lab)),r_unit)

heat_ybrks <- if (heat_gene_n <= 30) (heat_gene_lab) else (heat_ybrks)

if (max(heat_data$Rank) >= 18000) {
rank_scale <- max(heat_data$Rank)/18000
heat_data$Rank <- ceiling(heat_data$Rank / rank_scale)
heat_ybrks <- ceiling(heat_ybrks/rank_scale)
}

heat_n <- 0

for ( reps in c("All","Rep") ) {

  if (reps=="All") {
    heat_rep_data <- heat_data[heat_data$replicate == "All",]
  } else {
    heat_rep_data <- heat_data[heat_data$replicate != "All",]
  }

heatmap <- ggplot(
  heat_rep_data,
  aes(x=Position,y=Rank,fill=heat)
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
  limits=c(0.5,max(heat_rep_data$Rank)+0.5),
  breaks=heat_ybrks,
  expand=c(0,0)
  ) +
scale_x_continuous(
  limits=heat_xlim,
  breaks=heat_xbrks,
  expand=c(0,0)
  ) +
labs(
#  subtitle=paste("\n","\n",treat," vs ",control,"\n",protocol_title,"\n","n=",heat_gene_n,sep="")
#   subtitle=paste("n=",heat_gene_n,sep="")
    subtitle=ifelse(reps=="All",
      paste(
        "Mean of ", rep_paired, "s","\n",
        ifelse(common=="","",paste(common," ",sep="")),treat," vs ",control,"\n",
        protocol_title," ",prefix,"\n",
        "n=", heat_gene_n,
        sep=""
      ),
      paste(
        "Per ",rep_paired,
        sep=""
      )
    )
  ) +
xlab(
  paste(feature, "(bps)")
  ) +
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
    size=NA
    ),
  axis.title.y = element_text(size=9),
  axis.title.x = element_text(size=9),
  legend.title = element_text(size=9),
  plot.subtitle = element_text(size=9),
  legend.position = "right",
  strip.text=element_text(face="bold"),
  strip.background=element_rect(colour="transparent",fill="transparent",linewidth=0.1),
  axis.text.y = element_text(colour="black",size=9),
  axis.text.x = if (length(heat_xbrks)>2) (element_text(angle = 45, vjust = 1, hjust=1,colour="black")) else (element_text(colour="black"))
  )

if (length(unique(heat_data$Sense))<=1) {
  heatmap <- heatmap + theme(aspect.ratio=max(c((11/6),6/(8/length(unique(heat_rep_data$replicate))))))
} else {
  heatmap <- heatmap + theme(aspect.ratio=max(c((5/6),3/(8/length(unique(heat_rep_data$replicate)))))) 
}

if (length(unique(heat_data$Sense))>1) {
  if (reps == "All") { 
    all_heatmap <- heatmap + facet_wrap(vars(Sense),ncol=1,strip.position="right")
  } else {
    rep_heatmap <- heatmap + facet_grid(Sense~replicate)
  }
} else if (reps != "All") {
  rep_heatmap <- heatmap + facet_wrap(vars(replicate),nrow=1,strip.position="top")
} else {
  all_heatmap <- heatmap
}


}
#if (heat_n == 1) {
#  heat_legend <- as_ggplot(get_legend(heatmap))
#  heatmap <- heatmap +
#    theme(
#      legend.position="none"
#    )

heatmap <- ggarrange(
  plotlist=list(all_heatmap,rep_heatmap),
  nrow=2,
  ncol=1,
  heights=c(5,4),
 # labels=c(
 #   paste("Mean of ", rep_paired, "s",sep=""),
 #   paste("Per",rep_paired),
  vjust=1.5
)


assign(paste(r,"_heat",sep=""),heatmap) 

}


#heat_chunks <- split(heat_ls,ceiling(seq_along(heat_ls)/2))

#heatmap <- list()
#for ( c in heat_chunks) {

#heatmap_ls <- lapply(c,get)
#heatmap_ls <- append(heatmap_ls,heat_legend)

#heatmap_chunk <- ggarrange(plotlist=heatmap_ls, ncol=length(heatmap_ls), width=c(rep(2.5,length(heatmap_ls)-1),1), nrow=1)

#heatmap <- append(heatmap,heatmap_chunk)
#}


#heat_ls <- lapply(heat_ls,get)

#heatmap <- ggarrange(plotlist=heat_ls,ncol=2,nrow=length(heat_ls)/2,labels="AUTO")

heatmap <- lapply(heat_ls,get)

heatmap_title <- paste(
  "Heatmaps of changes ", difference, " of ", feature_i, " in ", experiment, ".",
  sep="")

heatmap_h <- ifelse(length(heat_config$bias)<=1,7.2,7.8)

heatmap_caption <- paste(
  "Heatmaps representing fold changes in ", difference, " of ", feature_i, " mapped with at least ", heat_min_reads, " reads on average in ", experiment, ".",
  sep="")



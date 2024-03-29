for (r in heat_config$Ranking) {
expr_heat <- expr_heat %>% arrange(!!sym(r)) %>% mutate(!!paste(r,"_Rank",sep='') := 1:n()) %>% ungroup
}


heat_pc = c(0,25,50,75,100)
heat_data <- heat_df[heat_df$featureID %in% expr_heat$featureID,]
heat_data$root_name <- expr_heat$root_name[match(heat_data$featureID,expr_heat$featureID)]

heat_data$covered <- ifelse(heat_data$coverage==0,0,1)

heat_data <- heat_data %>% arrange(Sense,covered,coverage) %>% group_by(Sense,covered,featureID) %>% mutate(feat_cov_rank = (1:n())/n() ) %>% ungroup

heat_data <- heat_data %>% arrange(Sense,covered,coverage) %>% group_by(Sense,covered) %>% mutate(heat_cov_rank = (1:n())/n() ) %>% ungroup

heat_data$feat_cov_rank <- heat_data$feat_cov_rank * heat_data$covered

heat_data$heat_cov_rank <- heat_data$heat_cov_rank * heat_data$covered

heat_data$heat <- ifelse(heat_data$coverage > 0 & heat_data$coverage >= min_heat_cov & heat_data$heat_cov_rank >= min_heat_cov_pc & heat_data$feat_cov_rank >= min_heat_cov_pc, heat_data$heat,0)

heat_ls <- paste(heat_ranks,"_heat",sep="")
heat_gene_n <- length(unique(heat_data$featureID))

heat_data$Sense <- factor(heat_data$Sense, levels=c("Sense","Antisense"))

heat_scale <- quantile(abs(heat_data$heat[heat_data$heat!=0]),heat_scale_pc,na.rm=T)
heat_data$heat <- heat_data$heat/heat_scale
heat_data$heat <- ifelse(heat_data$heat >= 1, 1, ifelse(heat_data$heat <= -1, -1, heat_data$heat))

heat_lfcbrks_i <- heat_lfcbrks/heat_scale

for (r in heat_config$Ranking) {

heat_unit <- heat_config$unit[heat_config$Ranking==r]
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

heatmap <- ggplot(
  heat_data,
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
  limits=c(0.5,max(heat_data$Rank)+0.5),
  breaks=heat_ybrks,
  expand=c(0,0)
  ) +
scale_x_continuous(
  limits=heat_xlim,
  breaks=heat_xbrks,
  expand=c(0,0)
  ) +
labs(
  subtitle=paste("n=",heat_gene_n,sep="")
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
  strip.text=element_text(face="bold"),
  strip.background=element_rect(colour="white",fill="white",size=0.1),
  axis.text.y = element_text(colour="black",size=9),
  axis.text.x = if (length(heat_xbrks)>2) (element_text(angle = 45, vjust = 1, hjust=1,colour="black")) else (element_text(colour="black"))
  )

if (length(unique(heat_data$Sense))==2) {

heatmap <- heatmap + facet_wrap(vars(Sense),ncol=1,strip.position="top")

}

assign(paste(r,"_heat",sep=""),heatmap) 

}

heat_ls <- lapply(heat_ls,get)

heatmap <- ggarrange(plotlist=heat_ls,ncol=2,nrow=length(heat_ls)/2,labels="AUTO")

heatmap_title <- paste(
  "Heatmaps of changes ", difference, " of ", feature_i, " in ", experiment, ".",
  sep="")

heatmap_h <- 7.2

heatmap_caption <- paste(
  "Heatmaps representing fold changes in ", difference, " of ", feature_i, " mapped with at least ", heat_min_reads, " reads on average in ", experiment, ".",
  sep="")



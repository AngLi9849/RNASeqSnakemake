cor_r <- cor(x=lfc_cor$x_lfc,y=lfc_cor$y_lfc,method="pearson",use="complete.obs")

cor_x_uplim <- quantile(lfc_cor$x_lfc,lfc_pc_lim/100,na.rm=T)
cor_x_lowlim <- quantile(lfc_cor$x_lfc,1-(lfc_pc_lim/100),na.rm=T)
cor_x_len <- abs(cor_x_uplim - cor_x_lowlim)

cor_y_uplim <- quantile(lfc_cor$y_lfc,lfc_pc_lim/100,na.rm=T)
cor_y_lowlim <- quantile(lfc_cor$y_lfc,1-(lfc_pc_lim/100),na.rm=T)
cor_y_len <- abs(cor_y_uplim - cor_y_lowlim)

cor_x_upqt <- quantile(lfc_cor$x_lfc,(1-cor_qt),na.rm=T)
cor_x_lowqt <- quantile(lfc_cor$x_lfc,cor_qt,na.rm=T)

cor_y_upqt <- quantile(lfc_cor$y_lfc,(1-cor_qt),na.rm=T)
cor_y_lowqt <- quantile(lfc_cor$y_lfc,cor_qt,na.rm=T)

cor_x_mid <- quantile(lfc_cor$x_lfc,0.5,na.rm=T)
cor_y_mid <- quantile(lfc_cor$y_lfc,0.5,na.rm=T)

lfc_cor$x_dist <- if (cor_x_len <=0) (0) else ((lfc_cor$x_lfc-cor_x_mid)/cor_x_len)
lfc_cor$y_dist <- if (cor_y_len <=0) (0) else ((lfc_cor$y_lfc-cor_y_mid)/cor_y_len)
lfc_cor$dist <-  lfc_cor$x_dist^2 + lfc_cor$x_dist^2
 
lfc_cor$cor <- ifelse( 
  ((lfc_cor$x_lfc >= cor_x_upqt) & (lfc_cor$y_lfc >= cor_y_upqt)) | ((lfc_cor$x_lfc <= cor_x_lowqt) & (lfc_cor$y_lfc <= cor_y_lowqt)),"positive",
  ifelse(
    ((lfc_cor$x_lfc <= cor_x_lowqt) & (lfc_cor$y_lfc >= cor_y_upqt)) | ((lfc_cor$x_lfc >= cor_x_upqt) & (lfc_cor$y_lfc <= cor_y_lowqt)), "negative"
    ,"insig"))

lfc_cor$position <- ifelse( lfc_cor$y_lfc >= cor_y_upqt, "top", ifelse( lfc_cor$y_lfc <= cor_y_lowqt, "bottom","middle"))

lfc_cor$Colour <- unlist(lapply(lfc_cor$cor,function(x) {get(paste(x,"_col",sep=""))}))


if (cor_r >= r_threshold) {
sig_cor <- "positive"
r_lab_pos <- 0
} else if ((0-cor_r) >= r_threshold) {
sig_cor <- "negative"
r_lab_pos <- 0.7
} else {
sig_cor <- "insig"
r_lab_pos <- 0
}

lfc_cor <- lfc_cor %>% arrange(cor,position,dist) %>% group_by(across(all_of(c("cor","position")))) %>% mutate(dist_rank=n():1) %>% ungroup()
lfc_cor$Label <- ifelse( (lfc_cor$cor != "insig") & (lfc_cor$dist_rank <= cor_lab_n) , lfc_cor$feature_name , NA) 

if (sig_cor_only) {
  lfc_cor$Colour <- ifelse(lfc_cor$cor==sig_cor, lfc_cor$Colour, insig_col)
  lfc_cor$Label <- ifelse(lfc_cor$cor==sig_cor, lfc_cor$Label, NA)
}

cor_xlim <- c(cor_x_lowlim, cor_x_uplim)
cor_ylim <- c(cor_y_lowlim, cor_y_uplim)
cor_xbrks <- as.numeric(signif(c(cor_x_mid - 0.9*(cor_x_mid-cor_x_lowlim), 0, cor_x_mid + 0.9*(cor_x_uplim-cor_x_mid)),2))
cor_ybrks <- as.numeric(signif(c(cor_y_mid - 0.9*(cor_y_mid-cor_y_lowlim), 0, cor_y_mid + 0.9*(cor_y_uplim-cor_y_mid)),2))

lfc_cor$x_lfc <- ifelse(lfc_cor$x_lfc > cor_x_uplim, cor_x_uplim, lfc_cor$x_lfc)
lfc_cor$x_lfc <- ifelse(lfc_cor$x_lfc < cor_x_lowlim, cor_x_lowlim, lfc_cor$x_lfc)
lfc_cor$y_lfc <- ifelse(lfc_cor$y_lfc > cor_y_uplim, cor_y_uplim, lfc_cor$y_lfc)
lfc_cor$y_lfc <- ifelse(lfc_cor$y_lfc < cor_y_lowlim, cor_y_lowlim, lfc_cor$y_lfc)

for (r in rankings) {
  lfc_cor <- lfc_cor %>% arrange(!!sym(r)) %>% mutate(!!paste(r,"_Rank",sep='') := 1:n()) %>% ungroup
}

correlation <- ggplot(data = lfc_cor, aes(
  x=x_lfc, 
  y=y_lfc, 
)) + 
  stat_cor(
    method="pearson",
    label.x.npc = r_lab_pos,
    label.y.npc= 0.999,
    size=3,
    colour="black") +
  geom_hline(
    yintercept=0,
    alpha=0.3,
    linetype=5
  ) +
  geom_vline(
    xintercept=0,
    alpha=0.3,
    linetype=5
  ) +
  labs(
    subtitle=paste(group_title,"\n",i_heading," genes","\nn=",gene_n,sep="")
  ) +
  scale_x_continuous(limits=cor_xlim, breaks = cor_xbrks) +
  scale_y_continuous(limits=cor_ylim, breaks = cor_ybrks) +
  xlab(paste(titles[[t]], "log2 FC")) +
  ylab(paste(main,"log2 FC")) +  
  theme(panel.background=element_rect(fill="White",colour="white"),
        strip.text=element_text(face="bold"),
        strip.background=element_rect(colour="white",fill="white",size=0.1),
        panel.border=element_rect(fill=NA,colour="black",size=0.7),
        legend.background=element_rect(fill="White"),
        legend.key=element_rect(colour="white",fill="White"),
        axis.line=element_line(colour="black",size=0.1),
        axis.line.x.top=element_line(colour="black",size=0.1),
        axis.line.y.right=element_line(colour="black",size=0.1),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9)
  )

correlate <- correlation + 
  geom_point(
    size=0.5,
    color=lfc_cor$Colour,
    shape = ifelse( 
      ( (lfc_cor$x_lfc <= cor_x_lowlim) | 
        (lfc_cor$x_lfc >= cor_x_uplim) | 
        ( lfc_cor$y_lfc <= cor_y_lowlim) | 
        (lfc_cor$y_lfc >= cor_y_uplim) ) , 17, 19)
  ) +
  geom_text_repel(
    mapping=aes(label=Label),
    size= 3.0,
    color=  "grey20",
    fontface= 1,
    segment.color= "grey20",
    segment.alpha= 0.5,
    bg.color= NA,
    bg.r= NA,
    max.overlaps=Inf,
    force=10,
    force_pull=5
  )

assign(paste("correlation_",cor_n,sep=""), correlate)

cor_n <- cor_n + 1

cor_rank_pc <- c(0,25,50,75,100)

for (r in rankings) {


r_rank <- paste(r,"_Rank",sep='')
cor_rank_brks <- unlist(lapply(cor_rank_pc, function(p) { quantile(as.numeric(unlist(lfc_cor[,r_rank])),p/100,na.rm=T) } ) )
names(cor_rank_brks) <-  unlist(lapply(cor_rank_pc, function(p) {signif(quantile(as.numeric(unlist(lfc_cor[,r])),p/100,na.rm=T),2) } ) )

# Aborted Attempt with bias on true values
#cor_bias_range <- cor_bias_max - cor_bias_min
#cor_bias_min <- signif(cor_bias_min - 0.05*cor_bias_range,2)
#cor_bias_max <- signif(cor_bias_max + 0.05*cor_bias_range,2)
#cor_bias_range <- cor_bias_max - cor_bias_min 
#cor_bias_lowqt <- signif(cor_bias_min + 0.25*cor_bias_range,2)
#cor_bias_mid <- signif(cor_bias_min + 0.5*cor_bias_range,2)
#cor_bias_upqt <- signif(cor_bias_min + 0.75*cor_bias_range,2)

cor_rank_lim <- c(1, nrow(lfc_cor))

cor_bias <- correlation + 
  geom_point(
    aes_string(colour=r_rank),
    size=1,
    alpha=0.3,
    shape = ifelse(
      ( (lfc_cor$x_lfc <= cor_x_lowlim) |
        (lfc_cor$x_lfc >= cor_x_uplim) |
        ( lfc_cor$y_lfc <= cor_y_lowlim) |
        (lfc_cor$y_lfc >= cor_y_uplim) ) , 17, 19)
  ) +
  geom_smooth(
    method = "lm",
    size=0.8,
    fill="black",
    colour="black",
    alpha=0.2) +
  scale_colour_gradientn(
    name=r,
    colours = cor_bias_colours,
    limits=cor_rank_lim,
    breaks=cor_rank_brks,
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
      direction = 'vertical'
    )
  )

assign(paste("correlation_",cor_n,sep=""), cor_bias)

cor_n <- cor_n + 1

}


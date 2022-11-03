cor_r <- cor(x=lfc_cor$x_lfc,y=lfc_cor$y_lfc,method="pearson",use="complete.obs")

cor_x_uplim <- quantile(lfc_cor$x_lfc,lfc_pc_lim/100,na.rm=T)
cor_x_lowlim <- quantile(lfc_cor$x_lfc,1-(lfc_pc_lim/100),na.rm=T)
cor_x_len <- abs(cor_x_uplim - cor_x_lowlim)

cor_y_uplim <- quantile(lfc_cor$y_lfc,lfc_pc_lim/100,na.rm=T)
cor_y_lowlim <- quantile(lfc_cor$y_lfc,1-(lfc_pc_lim/100),na.rm=T)
cor_y_len <- abs(cor_y_uplim - cor_y_lowlim)

cor_x_upqt <- quantile(lfc_cor$x_lfc,0.75,na.rm=T)
cor_x_lowqt <- quantile(lfc_cor$x_lfc,0.25,na.rm=T)

cor_y_upqt <- quantile(lfc_cor$y_lfc,0.75,na.rm=T)
cor_y_lowqt <- quantile(lfc_cor$y_lfc,0.25,na.rm=T)

cor_x_mid <- quantile(lfc_cor$x_lfc,0.5,na.rm=T)
cor_y_mid <- quantile(lfc_cor$y_lfc,0.5,na.rm=T)

lfc_cor$x_dist <- if (cor_x_len <=0) (0) else ((lfc_cor$x_lfc-cor_x_mid)/cor_x_len)
lfc_cor$y_dist <- if (cor_y_len <=0) (0) else ((lfc_cor$y_lfc-cor_y_mid)/cor_y_len)
lfc_cor$dist <-  lfc_cor$x_dist^2 + lfc_cor$x_dist^2
 

if (cor_r >= r_threshold) {

lfc_cor$position <- ifelse((lfc_cor$x_lfc >= cor_x_upqt) & (lfc_cor$y_lfc >= cor_y_upqt),"top",ifelse((lfc_cor$x_lfc <= cor_x_lowqt) & (lfc_cor$y_lfc <= cor_y_lowqt),"bottom","middle"))

cor_col <- pos_col
} else if ((0-cor_r) >= r_threshold) {

lfc_cor$position <- ifelse((lfc_cor$x_lfc <= cor_x_upqt) & (lfc_cor$y_lfc >= cor_y_upqt),"top",ifelse((lfc_cor$x_lfc >= cor_x_lowqt) & (lfc_cor$y_lfc <= cor_y_lowqt),"bottom","middle"))

cor_col <- neg_col
} else {

lfc_cor$position <- "middle"

cor_col <- insig_col
}

lfc_cor <- lfc_cor %>% arrange(position,dist) %>% group_by(position) %>% mutate(dist_rank=n():1) %>% ungroup()
lfc_cor$label <- ifelse( (lfc_cor$position != middle) & (lfc_cor$dist_rank <= cor_lab_n) , lfc_cor$feature_name , NA) 

lfc_cor$Colour <- ifelse(
  (lfc_cor$position != "middle"),
  cor_col,
  insig_col,
)

cor_xlim <- c(cor_x_lowlim, cor_x_uplim)
cor_ylim <- c(cor_y_lowlim, cor_y_uplim)
cor_xbrks <- signif(c(cor_x_mid - 0.9*(cor_x_mid-cor_x_lowlim), cor_x_mid, cor_x_mid + 0.9*(cor_x_uplim-cor_x_mid)),2)
cor_ybrks <- signif(c(cor_y_mid - 0.9*(cor_y_mid-cor_y_lowlim), cor_y_mid, cor_y_mid + 0.9*(cor_y_uplim-cor_y_mid)),2)



correlation <- ggplot(data = lfc_cor, aes(
  x=ifelse(x_lfc >= cor_x_uplim, 0.99*cor_x_len + cor_x_lowlim, ifelse(x_lfc <= cor_x_lowlim, cor_x_uplim - 0.99*cor_x_len, x_lfc)), 
  y=ifelse(y_lfc >= cor_y_uplim, 0.99*cor_y_len + cor_y_lowlim, ifelse(y_lfc <= cor_y_lowlim, cor_y_uplim - 0.99*cor_y_len, y_lfc)),
)) + 
  geom_point(
    size=0.5,
    color=lfc_cor$Colour,
    shape = ifelse( 
      ( (lfc_cor$x_lfc <= cor_x_lowlim) | 
        (lfc_cor$x_lfc >= cor_x_uplim) | 
        ( lfc_cor$y_lfc <= cor_y_lowlim) | 
        (lfc_cor$y_lfc >= cor_y_uplim) ) , 17, 19)
  ) +
  stat_cor(
    method="pearson",
    label.x.npc = 0,
    label.y.npc= 1,
    size=3,
    colour="black") +
  geom_hline(
    yintercept=cor_y_mid,
    alpha=0.3,
    linetype=5
  ) +
  geom_vline(
    xintercept=cor_x_mid,
    alpha=0.3,
    linetype=5
  ) +
  labs(
    subtitle=paste(i_heading,"\nn=",gene_n,sep="")
  ) +
  scale_x_continuous(limits=cor_xlim, breaks = cor_xbrks) +
  scale_y_continuous(limits=cor_ylim, breaks = cor_ybrks) +
  geom_text_repel(
    mapping=aes(label=feature_name),
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
  ) + 
  xlab(paste(titles[[t]])) +
  ylab(paste(main)) +  
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

assign(paste("correlation_",cor_n,sep=""), correlation)

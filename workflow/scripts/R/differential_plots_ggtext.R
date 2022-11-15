file_i <- gsub("_"," ",paste(a, "vs", control,  i, change,sep=" "))

dir_i <- paste(dir,"/",ifelse(i=="","All",i),sep="")
dir.create(dir_i)  

# Set variables for titles and legends
title <- gsub("_"," ",paste(a,"vs",control, toTitleCase(i),toTitleCase(change),sep=" "))
title
capt <- gsub("_"," ",paste(splice, tolower(prefix), counting, sep=" "))
capt
norm <- gsub("_"," ",paste(spikein, tolower(normaliser), sep=" "))
norm
  
max_log10P <- max(expr_i$log10P[expr_i$log10P < Inf],na.rm = T)
expr_i <-
  expr_i %>% dplyr::mutate(
    log10P = ifelse(
    log10P == Inf,
    max_log10P*1.1,
    log10P
  )
)

#Plot ma and volcano plots for each group
maplot <- ggplot(data = expr_i, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(
    size = 0.5,
    color = ifelse(expr_i$padj < 0.1, ifelse(expr_i$log2FoldChange < 0, "blue", "orange"), "grey20"),
    alpha = 0.5,
  ) +
  geom_point(
    data=subset(expr_i, gene_name %in% goi),
    color="red",
    size=1.5,
    alpha=1
  ) +
  scale_x_log10() +
  geom_hline(
    yintercept = 0,
    alpha=0.3,
    linetype=5
  ) +
  geom_text_repel(
    mapping=aes(label=ifelse((gene_name %in% goi), as.character(gene_name),ifelse((p_rank<=ma_n | lfc_rank <= ma_n) & padj < 0.1 ,as.character(gene_name),NA))),
    size=ifelse(expr_i$gene_name %in% goi, 4.0, 3.0),
    color=ifelse(expr_i$gene_name %in% goi, "white", "grey20"),
    fontface=ifelse(expr_i$gene_name %in% goi, 2, 1),
    segment.color=ifelse(expr_i$gene_name %in% goi, "black", "grey20"),
    segment.alpha=ifelse(expr_i$gene_name %in% goi, 1, 0.5),
    bg.color=ifelse(expr_i$gene_name %in% goi, "black", NA),
    bg.r=ifelse(expr_i$gene_name %in% goi, 0.2, NA),
    max.overlaps=Inf,
    force=10,
    force_pull=5
  ) +
  scale_y_continuous(
    limits=c(min(expr_i$log2FoldChange)*1.05,max(expr_i$log2FoldChange)*1.05)
  ) + 
  geom_hline(
    yintercept=0,
    alpha=0.3,
    linetype=5
  ) + 
  labs(
    caption=paste(
    "<br><span style=\"font-size:",as.character(title_size),"px;\"><strong>",
    title," MA plot.", 
    "</span></strong><br><span style=\"font-size:",as.character(legend_size),"px;\">", 
    capt," of ",length(expr_i$gene_name)," ", toTitleCase(i), " ", feature, "s normalised by ", norm, ". ", descript, 
    "</span><br>",
    sep="")
  ) +
  xlab(paste("Mean Count of",counting,sep=" ")) +
  ylab("log2 Fold Change") +
  theme(
    panel.background=element_rect(fill="White",colour="white"), 
    strip.text=element_text(face="bold"), 
    strip.background=element_rect(colour="white",fill="white",size=0.1), 
    panel.border=element_rect(fill=NA,colour="black",size=0.7), 
    legend.background=element_rect(fill="White"), 
    legend.key=element_rect(colour="white",fill="White"), 
    axis.line=element_line(colour="black",size=0.1), 
    axis.line.x.top=element_line(colour="black",size=0.1), 
    axis.line.y.right=element_line(colour="black",size=0.1),
    plot.caption=element_textbox(hjust=0,width=unit(1,"npc"),lineheight=1.5),
    plot.caption.position = "plot"
  )

volcplot <- ggplot(data = expr_i, aes(x = log2FoldChange, y = log10P)) +
  geom_point(
    size = 0.5,
    color = ifelse(expr_i$padj < 0.1, ifelse(expr_i$log2FoldChange < 0, "blue", "orange"), "grey20"),
    alpha = 0.5,
  ) +
  geom_point(
    data=subset(expr_i, gene_name %in% goi),
    color="red",
    size=1.5,
    alpha=1
  ) +
  geom_text_repel(
    mapping=aes(label=ifelse((gene_name %in% goi), as.character(gene_name),ifelse((p_rank<=volc_n | lfc_rank <= volc_n) & padj < 0.1 ,as.character(gene_name),NA))),
    size=ifelse(expr_i$gene_name %in% goi, 4.0, 3.0),
    color=ifelse(expr_i$gene_name %in% goi, "white", "grey20"),
    fontface=ifelse(expr_i$gene_name %in% goi, 2, 1),
    segment.color=ifelse(expr_i$gene_name %in% goi, "black", "grey20"),
    segment.alpha=ifelse(expr_i$gene_name %in% goi, 1, 0.5),
    bg.color=ifelse(expr_i$gene_name %in% goi, "black", NA),
    bg.r=ifelse(expr_i$gene_name %in% goi, 0.2, NA),
    max.overlaps=Inf,force=10,force_pull=5
  ) +
  scale_y_continuous(
    limits=c(0.0,max_log10P*1.1),expand=c(0.0,max_log10P/70)
  ) +
  scale_x_continuous(
    limits=c(-max(abs(expr_i$log2FoldChange)),max(abs(expr_i$log2FoldChange)))
  ) + 
  geom_vline(
    xintercept=0,
    alpha=0.3, 
    linetype=5
  ) + 
  labs(
    caption=paste(
    "<br><span style=\"font-size:",as.character(title_size),"px;\"><strong>",
    title," Volcano plot.",
    "</span></strong><br><span style=\"font-size:",as.character(legend_size),"px;\">",
    capt," of ",length(expr_i$gene_name)," ", toTitleCase(i), " ", feature,"s normalised by " , norm, ". ", descript,
    "</span><br>",
    sep="")
  ) +
  xlab(
    paste("log2 Fold Change of",change)
  ) +
  ylab("-log10 Adjusted P-Value") +
  theme(
    panel.background=element_rect(fill="White",colour="white"), 
    strip.text=element_text(face="bold"), 
    strip.background=element_rect(colour="white",fill="white",size=0.1), 
    panel.border=element_rect(fill=NA,colour="black",size=0.7), 
    legend.background=element_rect(fill="White"), 
    legend.key=element_rect(colour="white",fill="White"), 
    axis.line=element_line(colour="black",size=0.1), 
    axis.line.x.top=element_line(colour="black",size=0.1),  
    axis.line.y.right=element_line(colour="black",size=0.1),
    plot.caption=element_textbox(hjust=0,width=unit(1,"npc"),lineheight=1.5),
    plot.caption.position = "plot"
  )

ggsave(file=paste(file_i," MA Plot.pdf",sep=""), path=dir_i,plot=maplot,height=9,width=12,dpi=400)
ggsave(file=paste(file_i," MA Plot.png",sep=""), path=dir_i,plot=maplot,height=9,width=12,dpi=400)

ggsave(file=paste(file_i," Volcano Plot.pdf",sep=""), path=dir_i,plot=volcplot,height=9,width=12,dpi=400)
ggsave(file=paste(file_i," Volcano Plot.png",sep=""), path=dir_i,plot=volcplot,height=9,width=12,dpi=400)





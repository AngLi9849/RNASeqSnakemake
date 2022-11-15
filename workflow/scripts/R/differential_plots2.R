# Plot ma and volcano plots for all genes identified

ma_plot <- ggplot(data = expr, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(
    size = 0.5,
    color = ifelse(expr$padj < 0.1, ifelse(expr$log2FoldChange < 0, "blue", "orange"), "grey20"),
    alpha = 0.5,
  ) +
  geom_point(data=subset(expr, gene_name %in% goi),color="red",size=1.5,alpha=1) +
  scale_x_log10() + 
  geom_hline(
    yintercept = 0, 
    alpha=0.3, 
    linetype=5
  ) + 
  geom_text_repel(
    mapping=aes(label=ifelse((gene_name %in% goi), as.character(gene_name),ifelse((total_p_rank<=ma_n | total_lfc_rank <= ma_n) & padj < 0.1 ,as.character(gene_name),NA))),
    size=ifelse(expr$gene_name %in% goi, 4.0, 3.0),
    color=ifelse(expr$gene_name %in% goi, "white", "grey20"),
    fontface=ifelse(expr$gene_name %in% goi, 2, 1),
    segment.color=ifelse(expr$gene_name %in% goi, "black", "grey20"),
    segment.alpha=ifelse(expr$gene_name %in% goi, 1, 0.5),
    bg.color=ifelse(expr$gene_name %in% goi, "black", NA),
    bg.r=ifelse(expr$gene_name %in% goi, 0.2, NA),
    max.overlaps=Inf,force=10,force_pull=5
  ) +
  scale_y_continuous(limits=c(min(expr$log2FoldChange)*1.05,max(expr$log2FoldChange)*1.05)) + 
  labs(caption=paste("<br><span style=\"font-size:15px;\"><b>",title," MA plot ", "</b><br>", caption," of ",length(expr$gene_name)," genes.", descript, "</span><br>",sep="")) +
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

total_max_log10P <- max(expr$log10P[expr$log10P < Inf],na.rm = T)
expr_volc <-
  expr %>% dplyr::mutate(
    log10P = ifelse(
    log10P == Inf,
    total_max_log10P*1.1,
    log10P
  )
)

volc_plot <- ggplot(data = expr_volc, aes(x = log2FoldChange, y = log10P)) +
  geom_point(
    size = 0.5,
    color = ifelse(expr$padj < 0.1, ifelse(expr$log2FoldChange < 0, "blue", "orange"), "grey20"),
    alpha = 0.5,
  ) +
  geom_point(data=subset(expr, gene_name %in% goi),color="red",size=1.5,alpha=1) +
  geom_text_repel(
    mapping=aes(label=ifelse((gene_name %in% goi), as.character(gene_name),ifelse((total_p_rank<=volc_n | total_lfc_rank <= volc_n) & padj < 0.1 ,as.character(gene_name),NA))),
    size=ifelse(expr$gene_name %in% goi, 4.0, 3.0),
    color=ifelse(expr$gene_name %in% goi, "white", "grey20"),
    fontface=ifelse(expr$gene_name %in% goi, 2, 1),
    segment.color=ifelse(expr$gene_name %in% goi, "black", "grey20"),
    segment.alpha=ifelse(expr$gene_name %in% goi, 1, 0.5),
    bg.color=ifelse(expr$gene_name %in% goi, "black", NA),
    bg.r=ifelse(expr$gene_name %in% goi, 0.2, NA),
    max.overlaps=Inf,force=10,force_pull=5
  ) +
  scale_y_continuous(limits=c(0,total_max_log10P*1.1),expand=c(0,total_max_log10P/70)) +
  scale_x_continuous(limits=c(-max(abs(expr_volc$log2FoldChange)),max(abs(expr_volc$log2FoldChange)))) + 
  geom_vline(xintercept=0 ,alpha=0.3, linetype=5) + 
  ggtitle(paste(title,"Volcano plot",sep=" "), subtitle=paste(subtitle1," of ",length(expr_volc$gene_name)," genes",subtitle2, sep="")) +
  xlab(paste("log2 Fold Change of",counting)) +
  ylab("-log10 Adjusted P-Value") +
  theme(panel.background=element_rect(fill="White",colour="white"), strip.text=element_text(face="bold"), strip.background=element_rect(colour="white",fill="white",size=0.1), panel.border=element_rect(fill=NA,colour="black",size=0.7), legend.background=element_rect(fill="White"), legend.key=element_rect(colour="white",fill="White"), axis.line=element_line(colour="black",size=0.1), axis.line.x.top=element_line(colour="black",size=0.1), axis.line.y.right=element_line(colour="black",size=0.1))

write.table(expr,file=paste(dir,"/",title," Deseq.tab",sep=""),sep="\t",quote=F)

ggsave(file=paste(title,"MA Plot.pdf",sep=" "), path=dir,plot=ma_plot,height=9,width=12,dpi=400)
ggsave(file=paste(title,"MA Plot.png",sep=" "), path=dir,plot=ma_plot,height=9,width=12,dpi=400)

ggsave(file=paste(title,"Volcano Plot.png",sep=" "), path=dir,plot=volc_plot,height=9,width=12,dpi=400)
ggsave(file=paste(title,"Volcano Plot.pdf",sep=" "), path=dir,plot=volc_plot,height=9,width=12,dpi=400)

# Exclude genes not annotated with a biotype
expr <- expr[is.na(expr$biotype)==F,]
expr_bt <- expr[expr$biotype %in% biotypes,]

# Analyse genes grouped by biotype and whether it is monoexonic
for (i in unique(expr_bt$group)) {
ma_title_i <- gsub("_"," ",paste(a, "vs", control,  i, change, "MA Plot",sep=" "))
volc_title_i <- gsub("_"," ",paste(a, "vs", control,  i, change, "Volcano Plot",sep=" "))

dir_i <- paste(dir,"/",i,sep="")
dir.create(dir_i)  
  
expr_i <- expr_bt[expr_bt$group==i,]
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
    caption=paste("<br><span style=\"font-size:15px;\"><b>",title," MA plot ", "</b><br>", caption," of ",length(expr$gene_name)," genes.", descript, "</span><br>",sep="")
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
    caption=paste("<br><span style=\"font-size:15px;\"><b>",title," Volcano plot ", "</b><br>", caption," of ",length(expr$gene_name)," genes.", descript, "</span><br>",sep="")
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

ggsave(file=paste(ma_title_i,".pdf",sep=""), path=dir_i,plot=maplot,height=9,width=12,dpi=400)
ggsave(file=paste(ma_title_i,".png",sep=""), path=dir_i,plot=maplot,height=9,width=12,dpi=400)

ggsave(file=paste(volc_title_i,".pdf",sep=""), path=dir_i,plot=volcplot,height=9,width=12,dpi=400)
ggsave(file=paste(volc_title_i,".png",sep=""), path=dir_i,plot=volcplot,height=9,width=12,dpi=400)

}



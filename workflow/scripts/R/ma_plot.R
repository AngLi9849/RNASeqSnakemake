ma <- ggplot(data = expr_i, aes(x = baseMean, y = ifelse(abs(log2FoldChange) <= lfc_max, log2FoldChange, ifelse(log2FoldChange < 0, -0.98*lfc_max, 0.98*lfc_max)))) +
  geom_point(
    size = 0.5,
    color = expr_i$colour,
    alpha = 0.5,
    shape = ifelse(abs(expr_i$log2FoldChange) <= lfc_max, 19, 17)
  ) +
  scale_x_log10() +
  geom_hline(
    yintercept = 0,
    alpha=0.3,
    linetype=5
  ) +
  scale_y_continuous(
    limits=c(-lfc_max,lfc_max)
  ) +
  geom_hline(
    yintercept=0,
    alpha=0.3,
    linetype=5
  ) +
  xlab(toTitleCase(paste("Mean", counting,sep=" "))) +
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
    axis.title = element_text(size=9)
  )

ma_plot <- ma +
  geom_point(
    data=subset(expr_i, gene_name %in% goi),
    color="red",
    size=1.5,
    alpha=1
  ) +
  ylab(paste("log2 Fold Change of", feature_i, difference,sep= " ")) +
  xlab(toTitleCase(paste("Mean", feature_i, counting,sep=" "))) +  
  geom_text_repel(
    mapping=aes(label=ifelse((gene_name %in% goi), as.character(gene_name),ifelse((p_rank<=ma_n | lfc_rank <= ma_n) & padj < sig_p ,as.character(gene_name),NA))),
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
  ) 

ma <- ma +
  geom_point(
    data=subset(expr_i, gene_name %in% goi),
    color="red",
    size=1,
    alpha=1
  ) +
  geom_text_repel(
    mapping=aes(label=ifelse((gene_name %in% goi), as.character(gene_name),NA)),
    size=3.0,
    color="white",
    fontface=1,
    bg.color="black",
    bg.r=0.17,
    segment.color="black",
    segment.alpha=1,
    max.overlaps=Inf,
    force=10,
    force_pull=5
  )

ma_caption <- paste(
 "Log2 fold change in ", difference, " of each ", feature_i, " is plotted against their mean normalised ", counting, ". Dot colours represent significant (p < " ,sig_p,") increases (", up_col, ") and decreases (", down_col, ").",sep=""
)

ma_plot_title <- paste(title_i, "MA Plot.")
ma_plot_caption <- paste(
  ma_caption, " Names of top ", ma_n, " most significantly increased and decreased ", feature_i, " are labelled. ", ifelse(length(goi) > 0, "Genes of particular interest are highlighted with red dot and labelled in bold.", ""),sep=""
)

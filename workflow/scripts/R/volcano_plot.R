volcano <- ggplot(data = expr_i, aes(x = ifelse(abs(log2FoldChange) <= lfc_max, log2FoldChange, ifelse(log2FoldChange < 0, -0.98*lfc_max, 0.98*lfc_max)), y = log10P)) +
  geom_point(
    size = 0.5,
    color = expr_i$colour,
    alpha = 0.5,
    shape = ifelse(abs(expr_i$log2FoldChange) <= lfc_max, 19, 17)
  ) +
  scale_y_continuous(
    limits=c(0.0,max_log10P*1.1),expand=c(0.0,max_log10P/70)
  ) +
  scale_x_continuous(
    limits=c(-lfc_max,lfc_max)
  ) +
  geom_vline(
    xintercept=0,
    alpha=0.3,
    linetype=5
  ) +
  xlab(
    paste("log2 Fold Change")
  ) +
  ylab("-log10 P-Value") +
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

volcano_plot <- volcano +
  geom_point(
    data=subset(expr_i, gene_name %in% goi),
    color="red",
    size=1.5,
    alpha=1
  ) +
  xlab(paste("log2 Fold Change of",feature_i,difference)) +
  ylab("-log10 Adjusted P-Value") +
  geom_text_repel(
    mapping=aes(label=ifelse((gene_name %in% goi), as.character(gene_name),ifelse((p_rank<=volc_n | lfc_rank <= volc_n) & padj < sig_p ,as.character(gene_name),NA))),
    size=ifelse(expr_i$gene_name %in% goi, 4.0, 3.0),
    color=ifelse(expr_i$gene_name %in% goi, "white", "grey20"),
    fontface=ifelse(expr_i$gene_name %in% goi, 2, 1),
    segment.color=ifelse(expr_i$gene_name %in% goi, "black", "grey20"),
    segment.alpha=ifelse(expr_i$gene_name %in% goi, 1, 0.5),
    bg.color=ifelse(expr_i$gene_name %in% goi, "black", NA),
    bg.r=ifelse(expr_i$gene_name %in% goi, 0.2, NA),
    max.overlaps=Inf,force=10,force_pull=5
  )

volcano <- volcano + 
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

volcano_caption <- paste(
  "Significance (-log10 adjusted p-value) of changes in ", difference, " of each ", feature_i, " is plotted against the log2 fold change value. Dot colours represent significant (p < " ,sig_p,") increases (", up_col, ") and decreases (", down_col, ").",sep=""
)

volcano_plot_title <- paste(title_i, "Volcano Plot.")
volcano_plot_caption <- paste(
  volcano_caption, " Names of top ", volc_n, " most significantly increased and decreased are labelled. ", ifelse(length(goi) > 0, "Genes of particular interest are highlighted with red dot and labelled in bold.", ""),sep=""
)

 


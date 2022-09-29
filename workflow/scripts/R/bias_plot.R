expr_bias1 <- expr_i
expr_bias1$change <- "All"
expr_bias1$abslog2FoldChange <- expr_bias1$log2FoldChange
expr_bias1$density <- 1

expr_bias <- expr_i[abs(expr_i$padj) < sig_p,]

expr_bias$change <- paste("Significantly", expr_bias$change)
expr_bias$abslog2FoldChange <- abs(expr_bias$log2FoldChange)



FCmin=quantile(expr_bias$log2FoldChange,0.001,na.rm=TRUE)
FCmax=quantile(expr_bias$log2FoldChange,0.999,na.rm=TRUE)

expr_bias$density <- 0.01/expr_bias$padj
expr_bias <- rbind(expr_bias1,expr_bias)

GC_plot <- ggplot(data = expr_bias, aes(x=expr_bias$GC, y = expr_bias$abslog2FoldChange)) +
  geom_point(
    size=0.5,
    color=ifelse(expr_bias$log2FoldChange>0,up_col,down_col),
    alpha=expr_bias$density
  ) +
  geom_smooth(
    method = "lm",
    size=0.8,
    fill="black",
    colour="black",
    alpha=0.2) + 
  stat_cor(
    method="pearson",
    label.x.npc = 0, 
    label.y.npc= 1,
    size=3,
    colour="black") +
  geom_hline(
    yintercept=0,
    alpha=0.3,
    linetype=5
  ) +
  geom_vline(
    xintercept=0.5,
    alpha=0.3,
    linetype=5
  ) + 
  scale_x_continuous(limits=c(0.25,0.75),breaks=c(0,0.25,0.5,0.75,1),labels=scales::percent) + 
  facet_wrap(
    ~factor(change,levels=c("All","Significantly Increased","Significantly Decreased")), scales="free"
  ) +
  xlab(paste(feature_i,"GC Content (%)",sep=" ")) +
  ylab("log2 Fold Change") +
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

# Length to Log2foldChange with P-Value Color scale

Length_plot <- ggplot(data = expr_bias, aes(x=expr_bias$Length, y = expr_bias$abslog2FoldChange)) +
  geom_point(
    size=0.5,
    color=ifelse(expr_bias$log2FoldChange>0,up_col,down_col),
    alpha=expr_bias$density
  ) +
  geom_smooth(
    method = "lm",
    size=0.8,
    fill="black",
    colour="black",
    alpha=0.2) + 
  stat_cor(
    method="pearson",
    label.x.npc = 0, 
    label.y.npc= 1,
    size=3,
    colour="black") +
  scale_x_log10() + 
  facet_wrap(
    ~factor(change,levels=c("All","Significantly Increased","Significantly Decreased")), scales="free"
  ) +
  xlab(paste(feature_i, "Length (bps)",sep=" ")) +
  ylab("|log2 Fold Change|") +
  theme(panel.background=element_rect(fill="White",colour="white"), 
        strip.text=element_text(face="bold"), 
        strip.background=element_rect(colour="white",fill="white",size=0.1), 
        panel.border=element_rect(fill=NA,colour="black",size=0.7), 
        legend.background=element_rect(fill="White"), 
        legend.key=element_rect(colour="white",fill="White"), 
        axis.line=element_line(colour="black",size=0.1), 
        axis.line.x.top=element_line(colour="black",size=0.1), 
        axis.line.y.right=element_line(colour="black",size=0.1),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9)
  )

#rpkm bias plot

# Length to Log2foldChange with P-Value Color scale

rpkm_plot <- ggplot(data = expr_bias, aes(x=expr_bias$rpkm, y = expr_bias$abslog2FoldChange)) +
  geom_point(
    size=0.5,
    color=ifelse(expr_bias$log2FoldChange>0,up_col,down_col),
    alpha=expr_bias$density
  ) +
  geom_smooth(
    method = "lm",
    size=0.8,
    fill="black",
    colour="black",
    alpha=0.2) + 
  stat_cor(
    method="pearson",
    label.x.npc = 0, 
    label.y.npc= 1,
    size=3,
    colour="black") +
  scale_x_log10() + 
  facet_wrap(
    ~factor(change,levels=c("All","Significantly Increased","Significantly Decreased")), scales="free"
  ) +
  xlab(paste(feature_i, "Mean Expression Levels (RPKM)",sep=" ")) +
  ylab("|log2 Fold Change|") +
  theme(panel.background=element_rect(fill="White",colour="white"), 
        strip.text=element_text(face="bold"), 
        strip.background=element_rect(colour="white",fill="white",size=0.1), 
        panel.border=element_rect(fill=NA,colour="black",size=0.7), 
        legend.background=element_rect(fill="White"), 
        legend.key=element_rect(colour="white",fill="White"), 
        axis.line=element_line(colour="black",size=0.1), 
        axis.line.x.top=element_line(colour="black",size=0.1), 
        axis.line.y.right=element_line(colour="black",size=0.1),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9)
  )

bias_title <- paste("Potential Biases in ", title_i, ".",sep="")

if (nrow(expr_bias)==0) {
bias <- ""
bias_caption <- paste("There is not any significant changes in ", difference, " of ", feature_i, ".",sep="")
} else {
bias <- ggarrange(plotlist=list(GC_plot,Length_plot,rpkm_plot),ncol=1,nrow=3,labels="AUTO")
bias_caption <- paste(
  "Absolute values of significant (p < ", sig_p, ") log2 fold increases(", up_col, ") and decreases(", down_col, ") in ", difference, " of ", feature_i, " are plotted against their (A) GC content, (B) length and (C) mean expression levels in RPKM. A linear line of regression is shown with its range of dispersion (standard errors). R and p values of the regression analyses are indicated.", sep=""
)
}

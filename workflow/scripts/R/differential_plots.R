sink(log)
sink(log, type="message")

# Activate libraries for differential plotting
library(ggrepel)
library(ggtext)
library(stringr)
library(officer)
library(rvg)

# Import Plotting Parameters
descript <- as.character(snakemake@params[["descript"]])

ma_n <- as.numeric(snakemake@config[["differential_plots"]][["ma_gene_name_numbers"]])
volc_n <- as.numeric(snakemake@config[["differential_plots"]][["volcano_gene_name_numbers"]])

up_col <- as.character(snakemake@config[["differential_plots"]][["up_colour"]])
down_col <- as.character(snakemake@config[["differential_plots"]][["down_colour"]])
insig_col <- as.character(snakemake@config[["differential_plots"]][["insignificant_colour"]])

sig_p <- as.numeric(snakemake@config[["differential_plots"]][["significant_p"]])
uncert_p <- as.numeric(snakemake@config[["differential_plots"]][["uncertain_p"]]) 

plot_dpi <- as.numeric(snakemake@config[["differential_plots"]][["dpi"]])
tick_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["axis_ticks"]])
axis_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["axis_title"]])
label_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["labels"]])
title_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["title"]])
legend_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["legend"]])
font_size <- legend_size

ppt_w <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["width"]])
ppt_h <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["height"]])
ppt_top <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["margins"]][["top"]])
ppt_bottom <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["margins"]][["bottom"]])
ppt_left <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["margins"]][["left"]])
ppt_right <- as.numeric(snakemake@config[["differential_plots"]][["powerpoint"]][["margins"]][["right"]])

plot_w <- ppt_w - ppt_left - ppt_right
plot_h <- ppt_h - ppt_top - ppt_bottom

# Import word docx template and set fonts
doc <- read_docx(snakemake@input[["docx"]])

heading_1 <- fp_text(bold=TRUE,font.size=14)
heading_2 <- fp_text(bold=TRUE, font.size=13)
heading_3 <- fp_text(bold=TRUE,font.size=12)
plain <- fp_text(font.size=font_size)
bold <- fp_text(bold=TRUE, font.size=font_size)
italic <- fp_text(italic=TRUE, font.size=font_size)
italic_bold <- fp_text(bold=TRUE, italic=TRUE, font.size=font_size)

# Add heading to this differential analysis
heading <- gsub("_"," ",paste(a,"vs",control, feature, analysis ,sep=" "))
doc <- body_add(doc,fpar(ftext(heading,prop=heading_2)), style="heading 2")

# Plot figures for features in each mono/multiexonic-biotype groups
i_group <- append(unique(expr$group[expr$biotype %in% biotypes]), "")
for (i in i_group) {

if (i =="") {
  expr_i <- expr
} else {
  expr_i <- expr[expr$group==i,]
}

# Set file name and path
file_i <- gsub("_"," ",paste(a, "vs", control,  i, change,sep=" "))

dir_i <- paste(dir,"/",ifelse(i=="","All",i),sep="")
dir.create(dir_i)

# Set variables for titles and legends
title <- gsub("_"," ",paste(a,"vs",control, toTitleCase(i),toTitleCase(change),sep=" "))
capt <- gsub("_"," ",paste(splice, tolower(prefix), counting, sep=" "))
norm <- gsub("_"," ",paste(spikein, tolower(normaliser), sep=" "))

# Write heading for this analysis group
group_heading <- gsub("_"," ",toTitleCase(ifelse(i=="","Overview",i)))
doc <- body_add(doc,fpar(ftext(group_heading, prop=heading_3)),style = "heading 3")

# Restart figure counting
fig_num <- run_autonum(seq_id = "Figure", pre_label = "Figure ", post_label = ".", prop=bold,tnd=3, tns=".", bkm = "plot", start_at = 1)

plot_n <- 1

# Shrink Infinite log10Ps to 1.1 x maximum non-Inf log10P
max_log10P <- max(expr_i$log10P[expr_i$log10P < Inf],na.rm = T)
expr_i <-
  expr_i %>% dplyr::mutate(
    log10P = ifelse(
    log10P == Inf,
    max_log10P*1.1,
    log10P
  )
)

expr_i <- expr_i %>% arrange(change,padj) %>% group_by(change) %>% mutate(p_rank=1:n()) %>% ungroup
expr_i <- expr_i %>% arrange(change,abs(log2FoldChange)) %>% group_by(change) %>% mutate(lfc_rank=n():1) %>% ungroup

source(snakemake@config[["differential_plots"]][["scripts"]][["summary_plots"]])

plot_name <- "MA Plot"
expr_i$colour <- ifelse(expr_i$padj < p_threshold, ifelse(expr_i$log2FoldChange < 0, down_col, up_col), insig_col)
source(snakemake@config[["differential_plots"]][["scripts"]][["ma_plot"]])

plot_name <- "Volcano Plot"
expr_i$colour <- ifelse(expr_i$padj < p_threshold, ifelse(expr_i$log2FoldChange < 0, down_col, up_col), insig_col)
source(snakemake@config[["differential_plots"]][["scripts"]][["volcano_plot"]])

}

print(doc, target = snakemake@output[["docx"]])


sink(log)
sink(log, type="message")

# Import Differential Plotting Parameters
label_sum_plots <- as.logical(snakemake@config[["differential_plots"]][["label_summary"]])
ma_n <- as.numeric(snakemake@config[["differential_plots"]][["ma_gene_name_numbers"]])
volc_n <- as.numeric(snakemake@config[["differential_plots"]][["volcano_gene_name_numbers"]])
meta_trim <- as.numeric(snakemake@config[["metagene"]][["anomaly_trim"]]) 
lfc_pc_lim <- as.numeric(snakemake@config[["differential_plots"]][["lfc_pc_limit"]])

up_col <- as.character(snakemake@config[["differential_plots"]][["up_colour"]])
down_col <- as.character(snakemake@config[["differential_plots"]][["down_colour"]])
insig_col <- as.character(snakemake@config[["differential_plots"]][["insignificant_colour"]])
background_col <- "white"
void_col <- "black"

sig_p <- as.numeric(snakemake@config[["differential_analysis"]][["significant_p"]])
undetect_p <- as.numeric(snakemake@config[["differential_analysis"]][["undetect_p"]])

# Correlation settings
r_threshold <- snakemake@config[["group_analysis"]][["r_threshold"]]
cor_lab_n <- snakemake@config[["group_analysis"]][["label_number"]]
cor_page_n <- snakemake@config[["group_analysis"]][["correlation_per_page"]]
cor_qt <- as.numeric(snakemake@config[["group_analysis"]][["correlation_highlight_pc"]])/100
positive_col <- snakemake@config[["group_analysis"]][["positive_colour"]] 
negative_col <- snakemake@config[["group_analysis"]][["negative_colour"]]
sig_cor_only <- as.logical(snakemake@config[["group_analysis"]][["highlight_sig_cor_only"]])
cor_bias_colours <- lapply(strsplit(as.character(snakemake@config[["group_analysis"]][["cor_bias_colours"]]),","),trimws)[[1]]

# Heatmap setting
heat_name <- "log2FC"
heat_colours <- lapply(strsplit(as.character(snakemake@config[["heatmap"]][["heat_colours"]]),","),trimws)[[1]]
heat_lfcbrks_val <- c( -3 , -2 , -1 , -0.5 , 0 , 0.5 , 1 , 2 , 3 )
heat_lfcbrks <- unlist(lapply(heat_lfcbrks_val, function(x) {((2^(x+1))/((2^x)+1))-1}))
names(heat_lfcbrks) <- heat_lfcbrks_val
heat_scale_pc <- as.numeric(snakemake@config[["heatmap"]][["heat_scale_pc"]])/100

heat_ranks <- c("log2FoldChange","RPKM","Length","GC")
heat_units <- c("","","bps","%")

heat_config <- data.frame(heat_ranks,heat_units)
names(heat_config) <- c("Ranking","unit")
min_heat_cov <- snakemake@config[["heatmap"]][["min_cov"]]
min_heat_cov_pc <- snakemake@config[["heatmap"]][["min_cov_pc"]]

# General plotting settings
plot_dpi <- as.numeric(snakemake@config[["differential_plots"]][["dpi"]])
tick_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["axis_ticks"]])
axis_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["axis_title"]])
label_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["labels"]])
title_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["title"]])
legend_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["legend"]])
caption_size <- as.numeric(snakemake@config[["differential_plots"]][["font_sizes"]][["caption"]])
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
title_prop <- fp_text(font.size=title_size)
title_bold <- fp_text(bold=TRUE,font.size=title_size)
caption_prop <- fp_text(font.size=caption_size)
caption_bold <- fp_text(bold=TRUE,font.size=caption_size)

# A function to bold bracketed single capital letters for word output
format_captions <- function(c) {
  c_match <- str_match_all(c,"([^)]*\\()([^)]*)(\\)[^(]*)")
  c_list <- lapply(1:length(c_match), function(x) {
    if(nrow(c_match[[x]]) > 0) {
      unlist(apply(c_match[[x]], 1,function(y) {
        lapply(y[2:4],function(z){
          if(grepl("^[A-Z]$",z)[1]==TRUE){
            ftext(z,prop=caption_bold)
          } else {
            ftext(z,prop=caption_prop) 
          }
        })
      }),recursive = FALSE)
    } else {
      list(ftext(c[x],prop=caption_prop))
    }
  })
  return(c_list)
}


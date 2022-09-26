sink(log)
sink(log, type="message")

# Import Plotting Parameters

ma_n <- as.numeric(snakemake@config[["differential_plots"]][["ma_gene_name_numbers"]])
volc_n <- as.numeric(snakemake@config[["differential_plots"]][["volcano_gene_name_numbers"]])
meta_trim <- as.numeric(snakemake@config[["metagene"]][["anomaly_trim"]]) 

up_col <- as.character(snakemake@config[["differential_plots"]][["up_colour"]])
down_col <- as.character(snakemake@config[["differential_plots"]][["down_colour"]])
insig_col <- as.character(snakemake@config[["differential_plots"]][["insignificant_colour"]])
background_col <- "white"
void_col <- "black"

sig_p <- as.numeric(snakemake@config[["differential_analysis"]][["significant_p"]])
undetect_p <- as.numeric(snakemake@config[["differential_analysis"]][["undetect_p"]]) 

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


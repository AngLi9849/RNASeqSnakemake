#sig_bg <- read.csv(snakemake@input[["sig_bg"]],header=T,row.names = 1, sep='\t', check.names=FALSE)
sig <- as.numeric(snakemake@params[["sig"]])
bg <- as.numeric(snakemake@params[["bg"]])
section <- snakemake@params[["section"]]
base <- gsub("_"," ",gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@params[["base"]])))
base_feat <- gsub("_"," ",gsub("([^\\s_])([[:upper:]])([[:lower:]])",perl=TRUE,"\\1 \\2\\3",as.character(snakemake@params[["base_feat"]])))

mx_samples <- c(snakemake@params[["samples"]])

use_base_length <- (section!="body" & as.logical(snakemake@params[["main_int"]]))
use_base_length

# Import snakemake parameters and inputs
# Identify analysis
diff <- as.character(snakemake@wildcards[["difference"]])
difference <- as.character(diff)

analysis <- paste("differential", difference, "analysis")

biases <- read.csv(snakemake@config[["bias"]],header=T, sep='\t', check.names=FALSE)
biases
paste(biases$bias,"_count_bias",sep='')

if (difference != "splicing_ratio") {

  is_antisense <- snakemake@params[["is_antisense"]]==-1
  difference <- gsub("[_.]"," ",diff)
  difference_unit <- paste(difference, "(RPKM)")
  sig_bg <- read.csv(snakemake@input[["sig_bg"]],header=T,row.names = 1, sep='\t', check.names=FALSE)

  measurement <- "Count"

#  mx_df <- read.csv(snakemake@input[["mx_data"]],header=T, sep='\t', check.names=FALSE)
  bef_bin <- snakemake@params[["bef_bin"]]
  main_bin <- snakemake@params[["main_bin"]]
  section <- snakemake@params[["section"]]
  plot_median <- as.logical(snakemake@config[["metagene"]][["plot_median"]])
#  meta_y <- ifelse(plot_median,"Median","Mean")
  start_name <- snakemake@params[["start_name"]]
  start_name <- ifelse(start_name=="nan","Start",start_name)
  end_name <- snakemake@params[["end_name"]]
  end_name <- ifelse(end_name=="nan","End",end_name)
  if (is_antisense) {
    plotbef_len <- snakemake@params[["plotaft_len"]]
    plotaft_len <- snakemake@params[["plotbef_len"]]
    plotbef_bin <- snakemake@params[["plotaft_bin"]]
    plotaft_bin <- snakemake@params[["plotbef_bin"]]
    len_bef_n <- snakemake@params[["len_aft"]]
    len_aft_n <- snakemake@params[["len_bef"]]
  } else {
    plotbef_len <- snakemake@params[["plotbef_len"]]
    plotaft_len <- snakemake@params[["plotaft_len"]]
    plotbef_bin <- snakemake@params[["plotbef_bin"]]
    plotaft_bin <- snakemake@params[["plotaft_bin"]]
    len_bef_n <- snakemake@params[["len_bef"]]
    len_aft_n <- snakemake@params[["len_aft"]]
  }
  meta_bin <- plotbef_bin + plotaft_bin + main_bin

plotbef_brk_len <- signif(plotbef_len*1.5,1)/2
plotaft_brk_len <- signif(plotaft_len*1.5,1)/2
plotbef_brk_pos <- 0-bef_bin-(signif(plotbef_bin*1.5,1)/2)
plotaft_brk_pos <- main_bin-bef_bin + (signif(plotaft_bin*1.5,1)/2)
plotbef_brk <- paste(ifelse(is_antisense,"+","-"), as.character(plotbef_brk_len), sep="")
plotaft_brk <- paste(ifelse(is_antisense,"-","+"), as.character(plotaft_brk_len), sep="")

if (section=="body") {
meta_xbrks <- c(
 if(plotbef_len>0) (plotbef_brk_pos) else NULL,
 0,
 main_bin,
 if(plotaft_len>0) (plotaft_brk_pos) else NULL
)

names(meta_xbrks) <- c(
  if(plotbef_len>0) (plotbef_brk) else NULL,
  start_name,
  end_name,
  if(plotaft_len>0) (plotaft_brk) else NULL
)

} else {
start_name <- ifelse(start_name=="Start",paste(base_feat,toTitleCase(section)),start_name)

bef_brk_len <- signif(len_bef_n*1.5,1)/2
bef_brk <- paste(ifelse(is_antisense,"+","-"), as.character(bef_brk_len), sep="")
bef_brk_pos <- floor(signif(bef_bin*1.5,1)/2)

aft_brk_len <- signif(len_aft_n*1.5,1)/2
aft_brk <- paste(ifelse(is_antisense,"-","+"), as.character(aft_brk_len), sep="")
aft_brk_pos <- floor(signif((main_bin-bef_bin)*1.5,1)/2)

meta_xbrks <-  c(
  if(len_bef_n>0) (0-bef_brk_pos) else NULL,
  0,
  if(len_aft_n>0) (aft_brk_pos) else NULL
)

names(meta_xbrks) <- c(
  if(len_bef_n>0)(bef_brk) else NULL,
  start_name,
  if(len_aft_n>0) (aft_brk) else NULL
)

}

heat_df <- read_delim(snakemake@input[["heat_data"]],col_names=T,delim='\t')
heat_x_max <- max(heat_df$Position)
heat_x_min <- min(heat_df$Position)
heat_bin <- heat_x_max - heat_x_min
heat_xbrks <- meta_xbrks*heat_bin/meta_bin
heat_xlim <- c(heat_x_min-0.5,heat_x_max+0.5)

} else {
  difference <- gsub("_"," ",difference)
  difference_unit <- gsub("_"," ",difference)
  measurement <- diff
}


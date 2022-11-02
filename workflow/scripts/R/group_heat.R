heat_data <- bind_rows(lapply(lfc_i,function(x) {data.frame(x[,c("heat","title","rank")])})
names(heat_data) <- c("heat","title","rank")

heat_data$title <- factor(heat_data$title,levels=unlist(titles))

heat_scale <- quantile(abs(heat_data$heat),heat_scale_pc,na.rm=T)
heat_data$heat <- heat_data$heat/heat_scale
heat_data$heat <- ifelse(heat_data$heat >= 1, 1, ifelse(heat_data$heat <= -1, -1, heat_data$heat))

heat_lfcbrks_i <- heat_lfcbrks/heat_scale



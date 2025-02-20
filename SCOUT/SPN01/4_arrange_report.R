library(ggplot2)
library(patchwork)
library(dplyr)

page_rds <- list.files(path = "plot_rds/", pattern = "page",full.names = T)
pdf("SPN01_report.pdf", width = 8, height = 10)
lapply(page_rds, function(file) {
  data <- readRDS(file)
  print(data)	   
  # Assuming the data is plottable, modify accordingly to your data structure
  #plot <- ggplot2::ggplot(data, ggplot2::aes(x = <your_x_var>, y = <your_y_var>)) +
  #  ggplot2::geom_line() + # Modify depending on the plot type
  #  ggplot2::ggtitle(paste0("Plot for: ", basename(file)))
  #print(plot)
})
dev.off()


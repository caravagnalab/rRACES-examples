rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seq_res <- readRDS("data/seq_results.rds")
samples <- c('SPN01_Sample_1','SPN01_Sample_2','SPN01_Sample_3')
n_samples <- length(samples)

### VAFS ###
for (sample in samples) {
  pdf(paste0("tissue/chromosomes/", sample, '_report.pdf'), width = 8, height = 10)
  for (chr in unique(seq_res$chr)) {
    p1 <- plot_VAF_gw(seq_res, sample, chromosomes = chr)
    p2 <- plot_BAF_gw(seq_res, sample, chromosomes = chr)
    p3 <- plot_DR_gw(seq_res, sample, chromosomes = chr)
    
    p <- p1 / p2 / p3
    p %>% print()  
  }
  dev.off()  
}

### MARGINALS ###
marginals <- lapply(unique(seq_res$chr), function(c) {
  print(c)
  plot_marginals(seq_res, c, colour_by = 'classes')[[1]] + 
    ggtitle(paste0("Chr", c)) +
    scale_color_manual(values = c('driver'="black", 'passenger'="#F08080",'pre-neoplastic'= "#5F9EA0"))
})

p <- wrap_plots(marginals, ncol = 3) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave("tissue/chromosomes/vaf_marginals.pdf", dpi=300, width = 16, height = 20, plot = p)

setwd("~/GitHub/rRACES-examples/SPN04")
rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seq_res <- readRDS("data/sequencing.rds")
### VAFS ###
sample <- 'A'
for (sample in c('A', 'B')) {
  print(sample)
  pdf(paste0("tissue/chromosomes/", sample, '_report.pdf'), width = 8, height = 10)
  p1 <- plot_VAF_gw(seq_res, sample)
  p2 <- plot_BAF_gw(seq_res, sample)
  p3 <- plot_DR_gw(seq_res, sample)

  p <- p1 / p2 / p3
  p %>% print()
  dev.off()
}

### MARGINALS ###
seq_res$chr %>% unique()
marginals <- lapply(unique(seq_res$chr), function(c) {
  print(c)
  plot_marginals(seq_res, c, colour_by = 'classes')[[1]] +
    ggtitle(paste0("Chr", c)) +
    scale_color_manual(values = c('driver'="black", 'passenger'="#F08080",'pre-neoplastic'= "#5F9EA0"))
})

p <- wrap_plots(marginals, ncol = 3) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave("tissue/chromosomes/vaf_marginals.pdf", dpi=300, width = 16, height = 20, plot = p)

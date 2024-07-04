library(dplyr)
library(rRACES)
library(tidyverse)
library(patchwork)

#seq_results <- readRDS('/Users/lucreziavaleriani/Desktop/orfeo_LTS/races/SPN03/results/seq_80X.RDS')
seq_results <- readRDS('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/seq_80X.RDS')
samples <- c('Sample.A', 'Sample.B', 'Sample.C', 'Sample.D')

g_seq <- seq_results %>% filter(classes=="germinal")
for (sample in samples) {
  print(sample)
  
  p1 <- plot_VAF(g_seq, sample)
  p2 <- plot_BAF(g_seq, sample)
  p3 <- plot_DR(g_seq, sample)

  p <- p1 / p2 / p3
  ggsave(paste0("plots/SPN03_",sample, '_report.png'),plot = p, dpi=300, width = 8, height = 10)
}

pdf("plots/SPN03_vaf_report.pdf", width = 8, height = 10)
s_seq <- seq_results %>% filter(classes!="germinal")
for (c in unique(s_seq$chr)) {
    p_marg <- plot_VAF_marginals(s_seq, chromosomes = c, samples = samples, labels = s_seq["classes"]) 
    p_hist <- plot_VAF_histogram(s_seq, chromosomes = c, samples = samples, labels = s_seq["classes"], cuts = c(0.02, 1)) + xlim(0,1)
    p <- (wrap_plots(p_marg, ncol = 3) + p_hist) + 
      plot_layout(guides = 'collect', design = 'ABC\nDEF\nGGG\nGGG') + plot_annotation(title = (paste("Chromosome", c))) & theme(legend.position = 'bottom')
    print(p)
}
dev.off()

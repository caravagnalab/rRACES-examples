library(dplyr)
library(rRACES)
library(tidyverse)
library(patchwork)

seq_results <- readRDS("monday_meeting/WGD_post/seq_res_new.rds")
phylo_forest <- load_phylogenetic_forest("monday_meeting/WGD_post/phylo_forest.sff")

samples <- phylo_forest$get_samples_info()[["name"]] %>% sort()


#g_seq <- seq_results %>% filter(classes=="germinal")
#for (sample in samples) {
#  print(sample)
#  
#  p1 <- plot_VAF(g_seq, sample)
#  p2 <- plot_BAF(g_seq, sample)
#  p3 <- plot_DR(g_seq, sample)
#
#  p <- p1 / p2 / p3
#  ggsave(paste0("plots/SPN03_",sample, '_report.png'),plot = p, dpi=300, width = 8, height = 10)
#}

pdf("monday_meeting/SPN01_vaf_report.pdf", width = 8, height = 10)
s_seq <- seq_results %>% filter(classes!="germinal")
for (c in unique(s_seq$chr)) {
    p_marg <- plot_VAF_marginals(s_seq, chromosomes = c, samples = samples, labels = s_seq["classes"]) 
    p_hist <- plot_VAF_histogram(s_seq, chromosomes = c, samples = samples, labels = s_seq["classes"], cuts = c(0.01, 1)) + xlim(0,1)
    p <- (wrap_plots(p_marg, ncol = 3) + p_hist) + 
      plot_layout(guides = 'collect', design = 'ABC\n###\nGGG\nGGG') + plot_annotation(title = (paste("Chromosome", c))) & theme(legend.position = 'bottom')
    print(p)
}
dev.off()

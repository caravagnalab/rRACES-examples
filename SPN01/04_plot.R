rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seq_res <- readRDS("data/seq_results_80X_with_error_paired.rds")
### VAFS ###
samples <- c('SPN01_Sample_1','SPN01_Sample_2','SPN01_Sample_3')
g_seq <- seq_res %>% filter(classes=="germinal")
for (sample in samples) {
  print(sample)
  #pdf(paste0(sample, '_report.pdf'), width = 8, height = 10)
  #pdf(paste0("tissue/chromosomes/", sample, '_report.pdf'), width = 8, height = 10)
  p1 <- plot_VAF(g_seq, sample)
  p2 <- plot_BAF(g_seq, sample)
  p3 <- plot_DR(g_seq, sample)

  p <- p1 / p2 / p3
  ggsave(paste0("plots/",sample, '_report.png'),plot = p, dpi=300, width = 8, height = 10)
  #p %>% print()
  #dev.off()
}

### MARGINALS ###
pdf("plots/chromosome_vaf_marginals_report.pdf", width = 16, height = 5)
seq_res$chr %>% unique()
s_seq <- seq_res %>% filter(classes!="germinal")
for (c in unique(seq_res$chr)) {
    print(c)
    p_marg <- plot_VAF_marginals(s_seq, chromosomes = c, samples = samples, labels = s_seq["classes"])
    p_hist <- plot_VAF_histogram(s_seq, chromosomes = c, samples = samples, labels = s_seq["classes"], cuts = c(0.02, 1))
    p <- (wrap_plots(p_marg, ncol=4)+p_hist)+ plot_layout(guides = 'collect') & theme(legend.position = 'bottom') & ggtitle(paste("Chromosome", c))
    #p_marg <- lapply(p_marg, function(p) p + ggtitle(paste("Chromosome", c)))
    #p <- wrap_plots(list(p_marg,p_hist),ncol = 3, nrow=2) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
    print(p)
}

dev.off()
#marginals <- lapply(unique(seq_res$chr), function(c) {
#			    print(c)
#			    p_marg <-plot_VAF_marginals(s_seq, chromosomes=c, samples = samples, labels = s_seq["classes"])
#			    p_marg <- lapply(p_marg, function(p) p + ggtitle(paste("Chromosome", c)))
#			    p <- wrap_plots(p_marg, ncol = 3) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
#			    print(p)
#			    #ggsave(paste0("chr_",c,"_vaf_marginals.pdf"), dpi=300, width = 16, height = 8,plot=p)
#})


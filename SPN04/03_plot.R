library(dplyr)
library(rRACES)
library(tidyverse)
library(patchwork)

#seq_results <- readRDS('/Users/lucreziavaleriani/Desktop/orfeo_LTS/races/SPN03/results/seq_80X.RDS')
seq_results <- readRDS('data/seq_80X.RDS')
samples <- c('A', 'B')

gw_plots_baf <- lapply(samples, function(s){
  rRACES::plot_BAF(seq_results, sample = s, cuts = c(0, 1))
})
baf <- patchwork::wrap_plots(gw_plots_baf, nrow = 4)

gw_plots_dr <- lapply(samples, function(s){
  rRACES::plot_DR(seq_results, sample = s)
})
dr <- patchwork::wrap_plots(gw_plots_dr, nrow = 4)

gw_plots_vaf <- lapply(samples, function(s){
  rRACES::plot_VAF(seq_results, sample = s)
})
vaf <- patchwork::wrap_plots(gw_plots_vaf, nrow = 4)
gw_plot <- patchwork::wrap_plots(baf, dr, vaf)
ggsave(filename = 'plots/seq_plot.png', dpi = 300, plot = gw_plot,  width = 410, height = 297, units = "mm")

seq_results <- seq_results %>%
  select(!starts_with("normal")) %>%
  filter(causes!="germinal")

# Histogram
hist_class <- rRACES::plot_VAF_histogram(seq_results,
                                         cuts = c(0.05, 1),
                                         labels = seq_results["classes"])

hist_cause <- rRACES::plot_VAF_histogram(seq_results,
                                         cuts = c(0.05, 1),
                                         labels = seq_results["causes"])
hist <- hist_class + hist_cause + patchwork::plot_layout(nrow =2)
ggsave(filename = 'plots/histogram.png', plot = hist, dpi = 300,  width = 410, height = 297, units = "mm")


# Marginals
marg_class <- patchwork::wrap_plots(rRACES::plot_VAF_marginals(seq_results,
                                                               chromosome = '6',
                                                               labels = seq_results["classes"]), guides = 'collect') & theme_bw() + theme(legend.position = 'bottom')
marg_cause <-patchwork::wrap_plots(rRACES::plot_VAF_marginals(seq_results,
                                                              chromosome = '6',
                                                              labels = seq_results["causes"]), guides = 'collect') & theme_bw() + theme(legend.position = 'bottom')

marg <- patchwork::wrap_plots(marg_class, marg_cause, nrow = 2)
ggsave(filename = 'plots/marginal.png', plot = marg, dpi = 300,  width = 210, height = 297, units = "mm")

library(tidyverse)
library(patchwork)
source('../plotting/plotting_functions.R')

seq_results <- readRDS('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/seq_80X.RDS')
samples <- c('Sample.A', 'Sample.B', 'Sample.C', 'Sample.D')

gw_plots_baf <- lapply(samples, function(s){
  plot_BAF_gw(seq_results, sample = s, cuts = c(0, 1))
})
baf <- patchwork::wrap_plots(gw_plots_baf, nrow = 4)

gw_plots_dr <- lapply(samples, function(s){
  plot_DR_gw(seq_results, sample = s, cuts = c(0, 1))
})
dr <- patchwork::wrap_plots(gw_plots_dr, nrow = 4)

gw_plots_vaf <- lapply(samples, function(s){
  plot_VAF_gw(seq_results, sample = s, cuts = c(0.05, 1))
})
vaf <- patchwork::wrap_plots(gw_plots_vaf, nrow = 4)
gw_plot <- patchwork::wrap_plots(baf, dr, vaf)
ggsave(filename = './plots/gw_plot.png', dpi = 300, plot = gw_plot, width = 14, height = 15, units = 'in')


plot_histogram_vaf(seq_results, 
                   cuts = c(0.05, 1), 
                   colour_by = 'classes')
ggsave(filename = './plots/hist_class.png', dpi = 300, width = 16, height = 8, units = 'in')

plot_histogram_vaf(seq_results, 
                   cuts = c(0.05, 1), 
                   colour_by = 'causes')
ggsave(filename = './plots/hist_cause.png', dpi = 300, width = 16, height = 8, units = 'in')


patchwork::wrap_plots(plot_marginals(seq_results, 
                                     chromosome = '1', 
                                     colour_by = 'classes'), guides = 'collect') & theme_bw() + theme(legend.position = 'bottom')
ggsave(filename = './plots/marginal_class.png', dpi = 300, width = 10, height = 7, units = 'in')

patchwork::wrap_plots(plot_marginals(seq_results, 
                                     chromosome = '1', 
                                     colour_by = 'causes'), guides = 'collect') & theme_bw() + theme(legend.position = 'bottom')
ggsave(filename = './plots/marginal_cause.png', dpi = 300, width = 10, height = 7, units = 'in')



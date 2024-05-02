library(tidyverse)
library(patchwork)
source('../plotting/plotting_functions.R')

seq_results <- readRDS('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/seq_80X.RDS')
samples <- c('Sample.A', 'Sample.B', 'Sample.C', 'Sample.D')

gw_plots_baf <- lapply(samples, function(s){
  plot_BAF_gw(seq_results, sample = s, cuts = c(0, 1))
  })
patchwork::wrap_plots(gw_plots_baf)

gw_plots_dr <- lapply(samples, function(s){
  plot_DR_gw(seq_results, sample = s, cuts = c(0, 1))
})
patchwork::wrap_plots(gw_plots_dr)


gw_plots_vaf <- lapply(samples, function(s){
  plot_VAF_gw(seq_results, sample = s, cuts = c(0.05, 1))
})
patchwork::wrap_plots(gw_plots_vaf)

plot_histogram_vaf(seq_results, 
                   cuts = c(0.05, 1), 
                   colour_by = 'classes')

patchwork::wrap_plots(plot_marginals(seq_results, 
                                     chromosome = '1', 
                                     colour_by = 'classes'))

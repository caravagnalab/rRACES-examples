rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("../plotting/plotting_functions.R")

seq_rr <- readRDS("data/sequencing_homogeneous_growth.rds")

chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
for (chrom in chromosomes){
	marginal_plots <- plot_marginals(seq_rr, chromosome = chrom)
	marginals <- wrap_plots(marginal_plots,guides = "collect") & theme(legend.position = "bottom")
	ggsave(paste0("tissue_homo/chromosomes/",chrom, "_marginal.png"), dpi = 400, width = 8, height = 4, plot = marginals)
}

## Genome wide plots
#samples <- c('A', 'B', 'C')
#
#gw_plots_baf <- lapply(samples, function(s){
#  plot_BAF_gw(seq_rr, sample = s, cuts = c(0, 1))
#})
#baf <- patchwork::wrap_plots(gw_plots_baf, nrow = 4)
#
#gw_plots_dr <- lapply(samples, function(s){
#  plot_DR_gw(seq_rr, sample = s, cuts = c(0, 1))
#})
#dr <- patchwork::wrap_plots(gw_plots_dr, nrow = 4)
#
#gw_plots_vaf <- lapply(samples, function(s){
#  plot_VAF_gw(seq_rr, sample = s, cuts = c(0.05, 1))
#})
#vaf <- patchwork::wrap_plots(gw_plots_vaf, nrow = 4)
#gw_plot <- patchwork::wrap_plots(baf, dr, vaf)
#ggsave(filename = 'tissue_homo/gw_plot.png', dpi = 300, plot = gw_plot, width = 14, height = 15, units = 'in')

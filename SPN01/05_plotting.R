rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("../plotting/plotting_functions.R")

seq_rr <- readRDS("data/sequencing_homogeneous_growth.rds") %>%
  dplyr::rename(chr = chromosome)

plot_DR_gw(seq_rr, "A", chromosomes = c("1"))
plot_BAF_gw(seq_rr, "A")
plot_VAF_gw(seq_rr, sample = "A", chromosomes = c("1", "2", "6", "8"))
plot_histogram_vaf(seq_rr, cuts = c(.05, .098))

plot_list <- plot_marginals(seq_rr, chromosome = "6")

rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils.R")

seed <- 12345
set.seed(seed)

phylo_forest <- rRACES::load_phylogenetic_forest("data/phylo_forest.sff")

# Simulate sequencing ####
#chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
chromosomes <- c('22')

seq_results <- parallel::mclapply(chromosomes, function(c) {
	simulate_seq(phylo_forest, coverage = 80, chromosomes = c, write_SAM = FALSE)
}, mc.cores = parallel::detectCores()) %>% do.call("bind_rows", .)
#seq_results <- simulate_seq(phylo_forest, coverage = 80)

saveRDS(seq_results, "data/sequencing.rds")

rm(list=ls())
library(ProCESS)
library(tidyverse)

phylo_forest <- load_phylogenetic_forest("data/phylo_forest.sff")

cov <- 80
s <- simulate_seq(phylo_forest, coverage = cov, chromosomes = "1", write_SAM = FALSE)
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

seq_results <- parallel::mclapply(chromosomes, function(c) {
  simulate_seq(phylo_forest, coverage = cov, chromosomes = c, write_SAM = FALSE)
}, mc.cores = 8) %>% do.call("bind_rows", .)

saveRDS(seq_results, file = paste0('data/seq_', cov, 'X.RDS'))

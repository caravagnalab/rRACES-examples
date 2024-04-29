rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seed <- 12345
set.seed(seed)

phylo_forest <- rRACES::load_phylogenetic_forest("data/phylo_forest.sff")

chromosomes <- paste0(seq(1:22))
chromosomes <- c(chromosomes,"X", "Y")

seq_results <- parallel::mclapply(chromosomes, function(c) {
        simulate_seq(phylo_forest, coverage = 80, chromosomes = c, write_SAM = FALSE)
}, mc.cores = parallel::detectCores()) %>% do.call("bind_rows", .)

saveRDS(object =seq_results ,file = paste0("data/sequencing_homogeneous_growth.rds"))


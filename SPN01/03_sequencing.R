rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seed <- 12345
set.seed(seed)

phylo_forest <- load_phylogenetic_forest("data/phylo_forest.sff")

chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

# build a basic Illumina sequencer model in which errors occur
# at rate 4e-3 per base
basic_seq <- new(BasicIlluminaSequencer, 4e-3)

#seq_results <- parallel::mclapply(chromosomes, function(c) {
#	simulate_seq(phylo_forest, coverage = 80, chromosomes = c, write_SAM = FALSE)
#}, mc.cores = parallel::detectCores()) %>% do.call("bind_rows", .)

#saveRDS(object =seq_results ,file = paste0("data/sequencing_homogeneous_growth.rds"))
# let us simulate a 2.5x sequencing of the four sample
# on the error-less sequencer

chromosomes <- c("5","20")
seq_results <- parallel::mclapply(chromosomes, function(c) {
       simulate_seq(phylo_forest, sequencer = basic_seq,
		    coverage = 80, chromosomes = c, write_SAM = TRUE)
}, mc.cores = parallel::detectCores()) %>% do.call("bind_rows", .)

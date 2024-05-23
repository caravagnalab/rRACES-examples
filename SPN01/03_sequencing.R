rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seed <- 12345
set.seed(seed)

phylo_forest <- load_phylogenetic_forest("data/phylo_forest.sff")

curr_dir = getwd()
setwd("/orfeo/cephfs/scratch/cdslab/shared/races")


# Simulate sequencing ####
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

# seq_res <- lapply(chromosomes, function(x) {
# 	simulate_seq(phylo_forest, coverage = 80, chromosomes = x, write_SAM = TRUE, update_SAM = TRUE, rnd_seed = seed)
# })
# seq_res <- bind_rows(seq_res)

seq_results <- parallel::mclapply(chromosomes, function(c) {
	simulate_seq(phylo_forest, coverage = 80, chromosomes = c) 
		     #write_SAM = TRUE, 
		     #output_dir = "/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations",rnd_seed = seed)
}, mc.cores = parallel::detectCores()) 
seq_results_final<- do.call("bind_rows", seq_results)
#seq_results <- simulate_seq(phylo_forest, coverage = 80)

saveRDS(seq_results_final, paste0(curr_dir, "/data/seq_results.rds"))
print("sequencing ended")

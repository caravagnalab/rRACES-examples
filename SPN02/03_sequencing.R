library(rRACES)
library(dplyr)
library(ggplot2)

phylo_forest <- load_phylogenetic_forest("test/SPN02/phylo_forest.sff")
# forest <- load_samples_forest("test/SPN02/forest.sff")

# plot_forest(forest) %>% 
#     annotate_forest(phylo_forest, samples = T, MRCAs = T,
#                 exposures = T, drivers=T, add_driver_label = T)
# ggsave("test/SPN02/mutations_forest.pdf", width = 15)
curr_dir = getwd()
setwd("/orfeo/cephfs/scratch/cdslab/shared/races")

seed <- 12345
set.seed(seed)

# Simulate sequencing ####
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

# seq_res <- lapply(chromosomes, function(x) {
# 	simulate_seq(phylo_forest, coverage = 80, chromosomes = x, write_SAM = TRUE, update_SAM = TRUE, rnd_seed = seed)
# })
# seq_res <- bind_rows(seq_res)
basic_seq <- new(BasicIlluminaSequencer, 4e-3)

seq_results <- parallel::mclapply(chromosomes, function(c) {
	simulate_seq(phylo_forest, coverage = 80, chromosomes = c, 
	write_SAM = FALSE, rnd_seed = seed, 
	sequencer = basic_seq)
}, mc.cores = parallel::detectCores()) 
# %>% do.call("bind_rows", .)
#seq_results <- simulate_seq(phylo_forest, coverage = 80)

saveRDS(seq_results, paste0(curr_dir, "/test/SPN02/sequencing_wg_80x_witherr.rds"))
print("sequencing ended")

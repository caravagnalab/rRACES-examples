library(rRACES)
library(dplyr)
library(ggplot2)

phylo_forest <- load_phylogenetic_forest("test/SPN02/phylo_forest.sff")
# forest <- load_samples_forest("test/SPN02/forest.sff")

# plot_forest(forest) %>% 
#     annotate_forest(phylo_forest, samples = T, MRCAs = T,
#                 exposures = T, drivers=T, add_driver_label = T)
# ggsave("test/SPN02/mutations_forest.pdf", width = 15)

seed <- 12345
set.seed(seed)

# Simulate sequencing ####
# chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
chromosomes <- c("1", "2")
seq_results <- parallel::mclapply(chromosomes, function(c) {
	simulate_seq(phylo_forest, coverage = 80, chromosomes = c, write_SAM = TRUE, rnd_seed = seed)
}, mc.cores = parallel::detectCores()) %>% do.call("bind_rows", .)
#seq_results <- simulate_seq(phylo_forest, coverage = 80)

saveRDS(seq_results, "test/SPN02/sequencing_wg.rds")
print("sequencing ended")

library(rRACES)
library(dplyr)
library(ggplot2)

phylo_forest <- load_phylogenetic_forest("/orfeo/cephfs/scratch/area/vgazziero/CDSlab/SPN02/results/phylo_forest.sff")
# forest <- load_samples_forest("test/SPN02/forest.sff")

# plot_forest(forest) %>% 
#     annotate_forest(phylo_forest, samples = T, MRCAs = T,
#                 exposures = T, drivers=T, add_driver_label = T)
# ggsave("test/SPN02/mutations_forest.pdf", width = 15)
# curr_dir = getwd()
setwd("/orfeo/cephfs/scratch/cdslab/shared/mutation_engine/")

seed <- 12345
set.seed(seed)

# sam_folder = "/orfeo/cephfs/scratch/area/vgazziero/CDSlab/SPN02/results/purity_100"
# dir.create(sam_folder, recursive = TRUE)
# setwd(sam_folder)
# Simulate sequencing ####
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
# chromosomes <- c("1", "3", "6", "7", "22")
# seq_res <- lapply(chromosomes, function(x) {
# 	simulate_seq(phylo_forest, coverage = 80, chromosomes = x, write_SAM = TRUE, update_SAM = TRUE, rnd_seed = seed)
# })
# seq_res <- bind_rows(seq_res)
basic_seq <- BasicIlluminaSequencer(4e-3)

# print("starting seq and writing SAMs in desidered folder")

seq_results <- parallel::mclapply(chromosomes, function(c) {
	simulate_seq(phylo_forest, coverage = 80, chromosomes = c, 
		# output_dir = sam_folder, 
		write_SAM = FALSE, 
		sequencer = basic_seq, 
		purity = 1, 
		with_normal_sample = TRUE)
}, mc.cores = parallel::detectCores()) 
# %>% do.call("bind_rows", .)
#seq_results <- simulate_seq(phylo_forest, coverage = 80)
# list.files(sam_folder)

saveRDS(seq_results, "/orfeo/cephfs/scratch/area/vgazziero/CDSlab/SPN02/results/SPN02_seq_80x.rds")
print("sequencing ended")



curr_dir <- getwd()


# load phylogenetic forest
phylo_forest <- load_phylogenetic_forest( paste(curr_dir, "/data/phylo_forest.sff", sep = "") )

basic_seq <- BasicIlluminaSequencer(4e-3)

seq_results <- simulate_seq(phylo_forest, sequencer = basic_seq, coverage = 80, write_SAM = FALSE)

saveRDS(seq_results, paste0(curr_dir, "/data/sequencing_wg_80x_witherr.rds"))
print("sequencing ended")

#data <- seq_to_long(seq_results)
#BASE_FACTOR <- 1e6



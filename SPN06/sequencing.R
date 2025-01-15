
.libPaths("/orfeo/cephfs/scratch/cdslab/ahaghighi/Rlibs")

library(rRACES)
library(ggplot2)
library(dplyr)

seed <- 2024
set.seed(seed)

setwd("/u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06")

curr_dir <- getwd()

# load phylogenetic forest
phylo_forest <- load_phylogenetic_forest( paste(curr_dir, "/data/phylo_forest.sff", sep = "") )

chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
#chromosomes <- c("7","12","17","19","X","Y")


basic_seq <- BasicIlluminaSequencer(4e-3)

seq_results <- parallel::mclapply(
  chromosomes, 
  function(c) {
    simulate_seq(
      phylo_forest = phylo_forest, sequencer = basic_seq, chromosomes = c, coverage = 200, write_SAM = FALSE, purity = 1, with_normal_sample = TRUE
    )
  }, 
  mc.cores = parallel::detectCores()
)

#seq_results_final <- do.call("bind_rows", seq_results)

saveRDS(seq_results, file = paste(curr_dir, "/data/SPN06_seq_200X.rds", sep = "") )

print("sequencing done!")





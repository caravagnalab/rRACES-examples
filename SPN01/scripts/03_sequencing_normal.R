rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seed <- 12345
set.seed(seed)

phylo_forest <- load_phylogenetic_forest("data/phylo_forest.sff")

#curr_dir = getwd()
setwd("/orfeo/cephfs/scratch/cdslab/shared/races")


# Simulate sequencing ####
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
#no_error_seq <- ErrorlessIlluminaSequencer()
basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose
seq_results <- parallel::mclapply(chromosomes, function(c) {
       simulate_normal_seq(phylo_forest, chromosomes = c, coverage = 30,write_SAM = TRUE, read_size =150,
                    sequencer = basic_seq, insert_size_mean = 350, insert_size_stddev = 10, output_dir = paste0("/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/FINAL_DATA/sequencing_30X_basic_error_paired_100_1normal"),
       update_SAM =TRUE)
}, mc.cores = parallel::detectCores())

seq_results_final<- do.call("bind_rows", seq_results)
saveRDS(seq_results_final, "/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/data/seq_results_30X_basic_error_100_1normal.rds")
#####################################

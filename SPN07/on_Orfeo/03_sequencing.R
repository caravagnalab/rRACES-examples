#rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seed <- 12345
set.seed(seed)

phylo_forest <- load_phylogenetic_forest("/orfeo/LTS/CDSLab/LT_storage/antonelloa/my_home/rRACES-examples/SPN07/on_Orfeo/phyloforest.sff")

#curr_dir = getwd()
setwd("/orfeo/cephfs/scratch/cdslab/shared/races")


# Simulate sequencing ####
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
basic_seq <- BasicIlluminaSequencer(4e-3)
seq_results <- parallel::mclapply(chromosomes, function(c) {
  simulate_seq(phylo_forest, chromosomes = c, coverage = 80,write_SAM = FALSE,
               sequencer = basic_seq, insert_size = 150)
}, mc.cores = parallel::detectCores())

seq_results_final<- do.call("bind_rows", seq_results)
saveRDS(seq_results_final,"/orfeo/LTS/CDSLab/LT_storage/antonelloa/my_home/rRACES-examples/SPN07/on_Orfeo/seq_results_80X_with_error_paired.rds")
#print("sequencing ended")

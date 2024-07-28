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
	simulate_seq(phylo_forest, chromosomes = c, coverage = 100,write_SAM = TRUE,read_size =150, 
		     sequencer = basic_seq, insert_size_mean = 350, insert_size_stddev = 10, , output_dir = "/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/FINAL_DATA/sequencing_100X_basic_error_paired_350_3tumor_new1",
	update_SAM =TRUE, with_normal_sample =FALSE)
}, mc.cores = 8)

seq_results_final<- do.call("bind_rows", seq_results)
saveRDS(seq_results_final, "/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/data/seq_results_100X_basic_error_350_3tumor.rds")
#####################################

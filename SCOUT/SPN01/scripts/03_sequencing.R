rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
#print(args)
coverage = as.double(args[1])
purity = as.double(args[2])
output_dir = paste0("/orfeo/fast/cdslab/vgazziero/rRACES_test/SPN01/","spn01_",coverage,"X","_",purity,"p")
dir.create(paste0(output_dir,"/data"),recursive=TRUE)
dir.create(paste0(output_dir,"/SAM_files"),recursive=TRUE)
seed <- 12345
set.seed(seed)

phylo_forest <- load_phylogenetic_forest("data/phylo_forest_WGD.sff")

#curr_dir = getwd()
setwd("/orfeo/cephfs/scratch/cdslab/shared/races")


# Simulate sequencing ####
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
#no_error_seq <- ErrorlessIlluminaSequencer()
basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose



#seq_results <- parallel::mclapply(chromosomes, function(c) {
#	simulate_seq(phylo_forest, chromosomes = c, coverage = 20,write_SAM = TRUE,read_size =150, 
#		     sequencer = basic_seq, insert_size_mean = 350, insert_size_stddev = 10, , output_dir = "/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/CHECK_PURITY/sequencing_20X_basic_error_paired_350_3tumor_new1",
#	update_SAM =TRUE, with_normal_sample =FALSE)
#}, mc.cores = 8)


#chromosomes <- c("7","5")
seq_results <- parallel::mclapply(chromosomes, function(c) {
					          simulate_seq(phylo_forest, chromosomes = c, coverage = coverage, purity= purity, write_SAM = TRUE,read_size =150,
							                            sequencer = basic_seq, insert_size_mean = 350, insert_size_stddev = 10, , output_dir = paste0(output_dir,"/","SAM_files"), 
										            update_SAM =TRUE, with_normal_sample =FALSE)
}, mc.cores = 8)


seq_results_final<- do.call("bind_rows", seq_results)
saveRDS(seq_results_final, paste0(output_dir,"/data/seq_results_",coverage,"X_basic_error_",purity,"_p.rds"))
#####################################

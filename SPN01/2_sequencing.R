rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seed <- 12345
set.seed(seed)
dir <- getwd()
#chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
#chromosomes <- c("22")

chromosomes <- c("3")
phylo_forest <- load_phylogenetic_forest("monday_meeting/WGD_post/phylo_forest.sff")


#setwd("/orfeo/cephfs/scratch/cdslab/shared/races")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/new_simulation/sim6")
basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose



seq_results <- parallel::mclapply(chromosomes, function(c) {
					  	simulate_seq(phylo_forest, chromosomes = c, coverage = 100,write_SAM = FALSE,read_size =150, 
							       		     sequencer = basic_seq, insert_size_mean = 350, insert_size_stddev = 10, 
									     update_SAM =FALSE, with_normal_sample =TRUE,output_dir=paste0(dir,"data"))
}, mc.cores = 8)

seq_results_final<- do.call("bind_rows", seq_results)
saveRDS(seq_results_final, paste0(dir,"/seq_res_new_chr3.rds"))

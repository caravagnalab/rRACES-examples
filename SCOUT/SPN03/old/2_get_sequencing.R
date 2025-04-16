library(ProCESS)
library(tidyverse)

setwd('/orfeo/cephfs/scratch/cdslab/shared/races/')
phylo_forest <- load_phylogenetic_forest("/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/phylo_forest.sff")

cov <- 100
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
basic_seq <- BasicIlluminaSequencer(4e-3)


seq_results <- parallel::mclapply(chromosomes, function(c) {
        simulate_seq(phylo_forest, 
                     coverage = cov, 
                     chromosomes = c, 
                     write_SAM = FALSE,
                     read_size = 150,
                     sequencer = basic_seq,
                     insert_size_mean = 350, 
                     insert_size_stddev = 10,
                     output_dir = '/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/ProCESS_SAM',
                     update_SAM = TRUE, 
                     with_normal_sample = FALSE
                     )
}, mc.cores = 8)
seq_results_final<- do.call("bind_rows", seq_results)

saveRDS(object = seq_results_final, 
        file = paste0('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/LAST_seq_', cov, 'X.RDS'))

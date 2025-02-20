library(tidyverse)

cov_downsampled = read.table("/orfeo/cephfs/scratch/area/vgazziero/CDSlab/rRaces/downsampling/depth_test_downsampled", header = T)
chromosomes = c(paste0("chr", seq(1:22)), "chrX", "chrY")

colnames(cov_downsampled) = c("chr", "pos", "coverage")
cov_downsampled %>% 
    mutate(coverage = as.numeric(coverage)) %>% 
    ggplot2::ggplot(aes(coverage)) %>%
    geom_histogram()
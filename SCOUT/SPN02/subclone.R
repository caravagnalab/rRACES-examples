library(tidyverse)
library(rRACES)
library(CNAqc)

data = "/orfeo/cephfs/scratch/cdslab/shared/races/test_data_plotting/"

phylo_forest <- load_phylogenetic_forest(paste0(data, "phylo_forest.sff"))
sample_forest <- load_samples_forest(paste0(data, "samples_forest.sff"))
seq_res = readRDS(paste0(data, "seq_res_all.rds"))

bulk_allele_SA = phylo_forest$get_bulk_allelic_fragmentation("Sample_A")

seq_res_long = seq_to_long(seq_res)
mutations_sample_A = seq_res_long %>% 
    filter(sample_name == "Sample_A")

Sample_A_CNA_filtered = bulk_allele_SA %>% 
    filter(ratio > 0.05) %>% 
    rename(from = begin, to = end, Major = major, CCF = ratio)

Sample_A_CNA_filtered = Sample_A_CNA_filtered %>% 
    group_by(chr, from, to) %>% 
    group_split()

Sample_A_CNA_filtered = lapply(Sample_A_CNA_filtered, function(x) {
    x = x %>% 
        arrange(desc(CCF)) %>% 
        mutate(Major_2 = nth(Major, 2), minor_2 = nth(minor, 2))
    x = x[1,] 
    return(x)
}) %>% bind_rows()

mutations_sample_A_filtered = mutations_sample_A %>% 
    relocate(to, .after = from) %>% 
    filter(classes != "germinal")

sample_A_cnaqc = CNAqc::init(mutations = mutations_sample_A_filtered, 
    cna = Sample_A_CNA_filtered, 
    purity = 1, 
    sample = "sample_A", 
    ref = "GRCh38")


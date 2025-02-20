library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

forest <- load_samples_forest("/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/samples_forest.sff")

setwd('/orfeo/cephfs/scratch/cdslab/shared/races/')
m_engine <- MutationEngine(setup_code = "GRCh38",
                           tumour_type= 'CLLE')

SNV_rate = 1e-8
CNA_rate = 1e-11

# Clone 1 
# NOTCH1p2514*fs*4 
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = SNV_rate,
                                        CNA = CNA_rate),
                    drivers = list("NOTCH1 FY357Y", 
                                   CNA(chr = "13", 
                                       chr_pos = 39500001, 
                                       len = 1.5e6,  
                                       type = "D"))
)

# Clone 2 
# KRAS G12D
m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV = SNV_rate, 
                                        CNA = CNA_rate),
                    driver = list("KRAS G12D")
)

# Clone 3
# Unknown
m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = SNV_rate, 
                                        CNA = CNA_rate)
)

# Signatures
# SBS1, SBS5, IDSXX
m_engine$add_exposure(c(SBS5 = 0.3, SBS1 = 0.7, ID5 = 1))
m_engine

phylo_forest <- m_engine$place_mutations(forest, 
                                         num_of_preneoplatic_SNVs = 800,
                                         num_of_preneoplatic_indels = 200)

# check CN
# phylo_forest$get_sampled_cell_CNAs() %>% distinct(chr, begin, end)
phylo_forest$get_bulk_allelic_fragmentation('Sample A') %>% filter(major != minor, ratio > 0.9)
phylo_forest$get_bulk_allelic_fragmentation('Sample B') %>% filter(major != minor, ratio > 0.9)
phylo_forest$get_bulk_allelic_fragmentation('Sample C') %>% filter(major != minor, ratio > 0.9)

phylo_forest$save("/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/phylo_forest.sff")

annot_forest <- plot_forest(forest) %>%
    annotate_forest(phylo_forest, 
                    samples = T, 
                    MRCAs = T,
                    exposures = T, 
                    drivers = T, 
                    add_driver_label = T)

exp_timeline <- plot_exposure_timeline(phylo_forest)

labels <- get_relevant_branches(forest)
sticks <- plot_sticks(forest, labels)

pl <- annot_forest + sticks + exp_timeline + plot_layout(nrow = 3, design = 'A\nA\nB\nB\nC')
ggplot2::ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/plots/SPN03_mutations.png', plot = pl, width = 210, height = 297, units = "mm", dpi = 300)


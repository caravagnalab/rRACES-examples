rm(list = ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils.R")

seed <- 12345
set.seed(seed)

treatment_info <- readRDS("data/treatment_info.rds")
forest <- load_samples_forest("data/samples_forest.sff")

###### DRIVERS ######
setwd("/orfeo/cephfs/scratch/cdslab/shared/races/")

SNV_Clone1 = SNV(chr="2", chr_pos=209113113, alt="A")
CNA_Clone2 = CNA(type='A', chr='6', chr_pos=25100000, len=1e7 )
SNV_Clone3 = SNV(chr='1', chr_pos=115256530, alt='T')

mu_SNV <- 1e-8
mu_CNA <- 1e-11

###### MUTATION ENGINE ######
m_engine <- MutationEngine(setup_code = "GRCh38")

m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA),
                    drivers = list(SNV_Clone1))

m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA),
                    drivers = list(CNA_Clone2))

m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA),
                    drivers = list(SNV_Clone3))


###### SIGNATURE ######

m_engine$add_exposure(coefficients = c(SBS5 = 0.5, SBS1 = 0.5)) # ID1 is missing
m_engine$add_exposure(time = treatment_info$treatment_start, c(SBS5 = 0.4, SBS1 = 0.4, SBS25 = 0.2))
m_engine$add_exposure(time = treatment_info$treatment_end, coefficients = c(SBS5 = 0.5, SBS1 = 0.5))

###### PHYLO FOREST ######

phylo_forest <- m_engine$place_mutations(forest,
                                         num_of_preneoplatic_SNVs = 800,
                                         num_of_preneoplatic_indels = 200)
phylo_forest$save("/u/cdslab/gsantacatterina/scratch/ProCESS-examples/SPN04/data/phylo_forest.sff")


all_SNV <- phylo_forest$get_sampled_cell_mutations() %>% as_tibble()
all_SNV %>%
  group_by(cause) %>%
  summarise(nPos = n_distinct(chr_pos)) %>%
  print()

all_SNV %>%
  group_by(class) %>%
  summarise(nPos = n_distinct(chr_pos)) %>%
  print()

annot_forest <- plot_forest(forest) %>%
  annotate_forest(phylo_forest,
                  samples = T,
                  MRCAs = T,
                  exposures = T,
                  drivers=T,
                  add_driver_label = T)

exp_timeline <- plot_exposure_timeline(phylo_forest)

labels <- get_relevant_branches(forest)
sticks <- plot_sticks(forest, labels)

pl <- annot_forest + sticks + exp_timeline + plot_layout(nrow = 3, design = 'A\nA\nB\nB\nC')
pl <- annot_forest
ggplot2::ggsave('/u/cdslab/gsantacatterina/scratch/ProCESS-examples/SPN04/plots/mutations.png', plot = pl, width = 210, height = 297, units = "mm", dpi=300)
ggplot2::ggsave('/u/cdslab/gsantacatterina/scratch/ProCESS-examples/SPN04/plots/mutations.pdf', plot = pl, width = 210, height = 297, units = "mm", dpi=300)

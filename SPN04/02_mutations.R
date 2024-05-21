rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils.R")

seed <- 12345
set.seed(seed)

###### DRIVERS ######

SNV_Clone1 = SNV(chr="2", chr_pos=209113113, alt="A")
CNA_Clone2 = CNA(type='A', chr='6', chr_pos=25100000, len=1e7 )
SNV_Clone3 = SNV(chr='1', chr_pos=115256530, alt='T')

mu_SNV <- 2e-8
mu_CNA <- 1e-11

###### MUTATION ENGINE ######
m_engine <- build_mutation_engine(setup_code = "GRCh38", context_sampling = 20)

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
treatment_info <- readRDS("data/treatment_info.rds")
m_engine$add_exposure(coefficients = c(SBS5 = 0.5, SBS1 = 0.5)) # ID1 is missing
m_engine$add_exposure(time = treatment_info$treatment_start, c(SBS5 = 0.4, SBS1 = 0.4, SBS25 = 0.2))
m_engine$add_exposure(time = treatment_info$treatment_end, coefficients = c(SBS5 = 0.5, SBS1 = 0.5))

###### PHYLO FOREST ######
sampled_phylogeny <- load_samples_forest("data/samples_forest.sff")
phylo_forest <- m_engine$place_mutations(sampled_phylogeny, 1000)

all_SNV <- phylo_forest$get_sampled_cell_mutations() %>% as_tibble()
all_SNV %>%
  group_by(cause) %>%
  summarise(nPos = n_distinct(chr_pos)) %>%
  print()

all_SNV %>%
  group_by(class) %>%
  summarise(nPos = n_distinct(chr_pos)) %>%
  print()

tree_plot <- plot_forest(sampled_phylogeny)
tree_plot <- annotate_forest(tree_plot, forest = phylo_forest, exposures = TRUE, MRCAs = FALSE, samples = FALSE, facet_signatures = TRUE, drivers = TRUE, add_driver_label = FALSE)

phylo_forest$save("data/phylo_forest.sff")
ggsave("tissue/tree.pdf", dpi=300, width = 12, height = 12, plot = tree_plot)

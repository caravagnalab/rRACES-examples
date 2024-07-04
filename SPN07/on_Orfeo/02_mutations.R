library(rRACES)
library(dplyr)
library(ggplot2)

# sshfs Orfeo:/ dati_Orfeo 

#library(ggpubr)
# setwd("~/dati_Orfeo/orfeo/cephfs/home/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo")
#setwd("~/dati_Orfeo/orfeo/cephfs/scratch/cdslab/shared/races/")
setwd("/orfeo/cephfs/scratch/cdslab/shared/races/")

m_engine <- MutationEngine(setup_code = "GRCh38")

## 1. Add Drivers
m_engine$add_mutant(mutant_name = "1",
                    passenger_rates = c(SNV = 1e-9, CNA = 1e-11),
                    drivers = list(list('PTEN R130*', allele=0))
)
m_engine$add_mutant(mutant_name = "2",
                    passenger_rates = c(SNV = 1e-9, CNA = 1e-11),
                    drivers = list(
                      CNA("D", "10", 87862638 , 109292 ,allele = 1)
                    )
)
m_engine$add_mutant(mutant_name = "3",
                    passenger_rates = c(SNV = 1e-9, CNA = 1e-11),
                    drivers = list(
                      list('NF1 Q1174*', allele=0)
                    )
)
# Convert mutation in protein coordinates to genome coordinates : https://bibliome.ai/GRCh38/gene/ATRX
m_engine$add_mutant(mutant_name = "4",
                    passenger_rates = c(SNV = 1e-9, CNA = 1e-11),
                    drivers = list(
                      SNV('X', 2719, 'T','C') # ARTX R907*
                    )
)
m_engine$add_mutant(mutant_name = "5",
                    passenger_rates = c(SNV = 5e-9, CNA = 1e-11),
                    drivers = list(
                      SNV('2', 1082, 'A','G') # MSH6 c.1082G>A	p.R361H
                    )
)
m_engine$add_mutant(mutant_name = "6",
                    passenger_rates = c(SNV = 5e-9, CNA = 1e-11),
                    drivers = list(
                      list('TP53 R248W')
                    )
)

## 2. Add exposures
m_engine$add_exposure(c(SBS1 = .8, SBS5 = .2))
m_engine$add_exposure(c(ID1 = 1))
m_engine$add_exposure(time = 7.85,
                      c(SBS11= 1)) # Chemotherapy active from 7.85 to 10.85
m_engine$add_exposure(time = 10.85, c(ID1 = 1))
# Hypermutant signatures : SBS6, SBS14, SBS15, SBS20, SBS21, and SBS44.
m_engine$add_exposure(time = 10.85, c(SBS1 = .1, SBS5 = .1,SBS6=.4,
                                      SBS21=.2, SBS44=.2))


samples_forest <- load_samples_forest("/u/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo/forest_sampling.sff")
#samples_forest <- load_samples_forest('~/dati_Orfeo/orfeo/cephfs/home/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo/forest_sampling.sff')
phylo_forest <- m_engine$place_mutations(samples_forest, 1000)

# phylo_forest$get_sampled_cell_mutations() %>% head()
# phylo_forest$get_sampled_cell_CNAs() %>% head()
# phylo_forest$get_germline_mutations() %>% head()
phylo_forest$save("/u/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo/phyloforest.sff")

#phylo_forest <- load_phylogenetic_forest("phyloforest.sff")
#seq_results <- simulate_seq(phylo_forest, coverage = 50)

















rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
# Set directories
# dir.create(path = "data", recursive = TRUE)

# load the samples forest from "samples_forest.sff" and store it in `forest`
forest <- load_samples_forest("data/samples_forest.sff")
# building a mutation engine by using the "GRCh38" set-up configuration
setwd("/orfeo/cephfs/scratch/cdslab/shared/races/")

m_engine <- MutationEngine(setup_code = "GRCh38")

mu_SNV = 1e-8
mu_CNA = 1e-11
##112707518-112846239 
CNA_Clone2 = CNA(type = "D", "5",
                 chr_pos = 107707518, len = 1e7,allele = 0)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),
                    drivers = list(list("APC R1450*", allele = 1)))
m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),
                    drivers = list(CNA_Clone2))

m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),
                    drivers = list("KRAS G12V"))

m_engine$add_mutant(mutant_name = "Clone 4",
                    passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),
                    drivers = list("PIK3CA E545K"))


# Mutational signatures
m_engine$add_exposure(
  time = 0,
  coefficients = c(
    SBS1 = 0.40,
    SBS5 = 0.20,
    SBS10a = 0.15,
    SBS10b = 0.15,
    ID1 = 1,
    SBS15 = 0.1
  )
)

print("Mutation engine created")
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/data/phylo_forest.sff")
print("Mutations placed")


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
ggsave("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/plots/SPN01_mutations.png",plot = pl, width = 210, height = 297, units = "mm", dpi = 300)

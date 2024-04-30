rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
# Set directories
# dir.create(path = "data", recursive = TRUE)

# load the samples forest from "samples_forest.sff" and store it in `forest`
forest <- load_samples_forest("data/samples_forest_homogeneous_growth.sff")
# building a mutation engine by using the "GRCh38" set-up configuration
m_engine <- build_mutation_engine(setup_code = "GRCh38")
m_engine <- build_mutation_engine(setup_code = "GRCh38", context_sampling = 50)

## Drivers for the tumors
SNV_Clone1 = SNV("5", 112839942, "T",allele = 1) ## APC R1450*
CNA_Clone2 = CNA(type = "D", "5",
                 chr_pos = 67522671, len = 1e7,allele = 0)
SNV_Clone3 = SNV("12", 25245350, "A")
SNV_Clone4 = SNV("3",179218303, "A") # NC_000003.12:179218302:G:A

# Mutation rates (passengers)
mu_SNV = 1e-8
mu_CNA = 1e-11

m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = mu_SNV),
                    drivers = list(SNV_Clone1))

m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(CNA = mu_CNA),
                    drivers = list(CNA_Clone2))

m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = mu_SNV),
                    drivers = list(SNV_Clone3))

m_engine$add_mutant(mutant_name = "Clone 4",
                    passenger_rates = c(SNV = mu_SNV),
                    drivers = list(SNV_Clone4))

# Mutational signatures
m_engine$add_exposure(
  time = 0,
  coefficients = c(
    SBS1 = 0.4,
    SBS5 = 0.3,
    SBS13 = 0.3)
)

phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_mutations = 1000)
phylo_forest$save("data/phylo_forest.sff")

chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

seq_results <- parallel::mclapply(chromosomes, function(c) {
	simulate_seq(phylo_forest, coverage = 80, chromosomes = c, write_SAM = FALSE)
}, mc.cores = parallel::detectCores()) %>% do.call("bind_rows", .)

saveRDS(object =seq_results ,file = paste0("data/sequencing_homogeneous_growth.rds"))


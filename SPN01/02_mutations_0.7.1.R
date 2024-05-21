rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
# Set directories
# dir.create(path = "data", recursive = TRUE)

# load the samples forest from "samples_forest.sff" and store it in `forest`
forest <- load_samples_forest("samples_forest_final.sff")
# building a mutation engine by using the "GRCh38" set-up configuration
m_engine <- build_mutation_engine(setup_code = "GRCh38")
m_engine <- build_mutation_engine(setup_code = "GRCh38", context_sampling = 50)

mu_SNV = 1e-8
mu_CNA = 1e-11
## APC gene coordinates 112707518-112846239 
## Delited region=107707518-117707518
CNA_Clone2 = CNA(type = "D", "5",
                 chr_pos = 107707518, len = 1e7,allele = 0)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = mu_SNV),
                    drivers = list(list("APC R1450*", allele = 1)))
m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(CNA = mu_CNA),
                    drivers = list(CNA_Clone2))

m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = mu_SNV),
                    drivers = list("KRAS G12V"))

m_engine$add_mutant(mutant_name = "Clone 4",
                    passenger_rates = c(SNV = mu_SNV),
                    drivers = list("PIK3CA E545K"))


# Mutational signatures
m_engine$add_exposure(
  time = 0,
  coefficients = c(
    SBS1 = 0.4,
    SBS5 = 0.3,
    SBS13 = 0.3)
)

print("Mutation engine created")
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_mutations = 1000)
phylo_forest$save("phylo_forest_final.sff")
print("Mutations placed")

chromosomes <- c("5")
seq_results <-simulate_seq(phylo_forest,chromosomes = c("5"),coverage = 40)
#seq_results <- parallel::mclapply(chromosomes, function(c) {
#       simulate_seq(phylo_forest,chromosomes = c,
#		    coverage = 80) 
#		    output_dir = "rRACES_SAM_1",update_SAM = TRUE,write_SAM = TRUE)
#}, mc.cores = parallel::detectCores())
##%>% do.call("bind_rows", .)
saveRDS(object =seq_results, file = "seq_results_5.rds")
seq_results_final<- do.call("bind_rows", seq_results)
saveRDS(object =seq_results_final, file = "seq_results.rds")

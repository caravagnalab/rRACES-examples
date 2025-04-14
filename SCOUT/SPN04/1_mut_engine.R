rm(list = ls())
library(ProCESS)
library(dplyr)
seed <- 12345
set.seed(seed)

outdir <- "/orfeo/scratch/cdslab/shared/SCOUT/SPN04/races/"
forest <- load_samples_forest(paste0(outdir,"sample_forest.sff"))
treatment_info <- readRDS("treatment_info.rds")

setwd("/orfeo/scratch/cdslab/shared/races/GRCh38/")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_study = "US")

mu_SNV <- 1e-8
mu_CNA <- 1e-11

SNV_Clone1 = SNV(chr="2", chr_pos=209113113, alt="A")
CNA_Clone2 = CNA(type='A', chr='6', chr_pos=25100000, len=1e7 )
SNV_Clone3 = SNV(chr='1', chr_pos=115256530, alt='T')



# Add drivers ####
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA),
                    drivers = list(SNV_Clone1))

m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA),
                    drivers = list(CNA_Clone2))

m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA),
                    drivers = list(SNV_Clone3))


# Mutational Signatures ####
m_engine$add_exposure(coefficients = c(SBS5 = 0.5, SBS1 = 0.5)) # ID1 is missing
m_engine$add_exposure(time = treatment_info$treatment_start, c(SBS5 = 0.4, SBS1 = 0.4, SBS25 = 0.2))
m_engine$add_exposure(time = treatment_info$treatment_end, coefficients = c(SBS5 = 0.5, SBS1 = 0.5))

# Phylo forest ######
print("Mutation engine created")
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save(paste0(outdir,"phylo_forest.sff"))
dir.create(paste0(outdir,"cna_data_v4"))
sample_names <- phylo_forest$get_samples_info()[["name"]]
lapply(sample_names,function(s){
  cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
  saveRDS(file=paste0(outdir,"cna_data_v4/",s,"_cna.rds"),object=cna)
})

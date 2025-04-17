rm(list = ls())
library(ProCESS)
library(dplyr)
seed <- 12345
set.seed(seed)

outdir <- "/orfeo/scratch/cdslab/shared/SCOUT/SPN04/process/"
forest <- load_samples_forest(paste0(outdir,"sample_forest.sff"))
treatment_info <- readRDS(paste0(outdir,"treatment_info.rds"))

setwd("/orfeo/scratch/cdslab/shared/ProCESS/GRCh38/")
# Set tumour type and gender
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_study = "KR", tumour_type = "LAML")
m_engine$set_germline_subject("NA18940")

mu_SNV = 1e-8
mu_ID = 1e-9
mu_CNA = 1e-11

driver_Clone1 = SNV(chr="2", chr_pos=208248389, alt="T", ref="C", allele=1)
driver_Clone2 = CNA(type='A', chr='12', chr_pos=25100000, len=1e7)
driver_Clone3 = list("NRAS Q61K", allele = 1)

# Add drivers ####

m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA, indel=mu_ID),
                    drivers = list(driver_Clone1))

m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA, indel=mu_ID),
                    drivers = list(driver_Clone2))

m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA, indel=mu_ID),
                    drivers = list(driver_Clone3))


# Mutational Signatures ####
m_engine$add_exposure(coefficients = c(SBS5 = 0.5, SBS1 = 0.5, ID1 = 1))
m_engine$add_exposure(time = treatment_info$treatment_start, c(SBS5 = 0.35, SBS1 = 0.35, SBS25 = 0.3, ID1 = 1))
m_engine$add_exposure(time = treatment_info$treatment_end, coefficients = c(SBS5 = 0.5, SBS1 = 0.5, ID1 = 1))

# Phylo forest ######
print("Mutation engine created")
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save(paste0(outdir,"phylo_forest.sff"))

dir.create(paste0(outdir,"cna_data"))
sample_names <- phylo_forest$get_samples_info()[["name"]]
lapply(sample_names,function(s){
  cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
  saveRDS(file=paste0(outdir,"cna_data/",s,"_cna.rds"),object=cna)
})

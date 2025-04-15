rm(list = ls())
library(ProCESS)
library(dplyr)
set.seed(06117)


outdir <- "/orfeo/scratch/cdslab/shared/SCOUT/SPN01/races/"
forest <- load_samples_forest(paste0(outdir,"sample_forest.sff"))


setwd("/orfeo/cephfs/scratch/cdslab/shared/ProCESS/GRCh38")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COAD",
                           tumour_study = "US")


mu_SNV = 1e-8
mu_CNA = 1e-10
##112707518-112846239 
CNA_Clone2 = ProCESS::CNA(type = "D", "5",
                         chr_pos = 107707518, len = 2e7,allele = 0)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = mu_SNV, CNA = 0),drivers = list(list("APC R1450*", allele = 1)))
m_engine$add_mutant(mutant_name = "Clone 2",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),drivers = list(CNA_Clone2))
mu_SNV = 1e-8
mu_CNA = 1e-13
m_engine$add_mutant(mutant_name = "Clone 3",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),drivers = list("TP53 R175H"))

m_engine$add_mutant(mutant_name = "Clone 4",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),drivers = list(WGD))


# Mutational signatures
m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.20,SBS5 = 0.35,
                                                SBS18 = 0.15,SBS17b = 0.25,ID1 = 0.60,ID2 = 0.40,SBS88 = 0.05))
print("Mutation engine created")
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save(paste0(outdir,"phylo_forest.sff"))
dir.create(paste0(outdir,"cna_data_v4"))
sample_names <- phylo_forest$get_samples_info()[["name"]]
lapply(sample_names,function(s){
    cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
    saveRDS(file=paste0(outdir,"cna_data_v4/",s,"_cna.rds"),object=cna)
})

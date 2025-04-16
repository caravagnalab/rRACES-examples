library(ProCESS)
library(dplyr)

base="/orfeo/cephfs/scratch/area/lvaleriani/races/SPN03_rRACES_report"
setwd(base)

forest <- load_samples_forest("sample_forest.sff")

setwd('/orfeo/cephfs/scratch/cdslab/shared/ProCESS/GRCh38')
m_engine <- MutationEngine(setup_code = "GRCh38",
                           tumour_type= 'CLLE')

SNV_rate = 1e-8
indel_rate = 1e-9
CNA_rate = 1e-12

# Clone 1 
# NOTCH1p2514*fs*4 
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = SNV_rate,
                                        CNA = CNA_rate,
                                        indel = indel_rate),
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
                                        CNA = CNA_rate,
                                        indel = indel_rate),
                    driver = list("KRAS G12D")
)

# Clone 3
# Unknown
m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = SNV_rate,
                                        CNA = CNA_rate,
                                        indel = indel_rate)
)

# Signatures
# SBS1, SBS5, IDSXX
m_engine$add_exposure(c(SBS5 = 0.3, SBS1 = 0.7, ID5 = 1))
m_engine

phylo_forest <- m_engine$place_mutations(forest, 
                                         num_of_preneoplatic_SNVs = 800,
                                         num_of_preneoplatic_indels = 200)

phylo_forest$save(paste0(base,"phylo_forest.sff"))


dir.create(paste0(base, "cna_data"),recursive = T)
sample_names <- phylo_forest$get_samples_info()[["name"]]
lapply(sample_names,function(s){
  cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
  saveRDS(file=paste0(base, "cna_data/",s,"_cna.rds"), object=cna)
})




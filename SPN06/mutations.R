
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

#-------------------------------------------------------------------------------
#-------------------------- set up Mutation Engine -----------------------------
#-------------------------------------------------------------------------------
curr_dir <- getwd()

#setwd("/orfeo/cephfs/scratch/cdslab/shared/mutation_engine/")

m_engine <- MutationEngine(setup_code = "GRCh38", tumour_type = "LUAD", tumour_study = "US")#, context_sampling = 50)

message("start creating mutations")
#-------------------------------------------------------------------------------
#---------------------------- passenger mutations ------------------------------
#-------------------------------------------------------------------------------
passengers <- c(SNV = 1e-8, CNA = 1e-11)
#SNV(chr = 0, chr_pos = 0, alt = 0, ref = "?", allele = NULL, cause = "")

#-------------------------------------------------------------------------------
#------------------------------ driver mutations -------------------------------
#-------------------------------------------------------------------------------

# TP53 R248Q/W/L
SNV_C1 = "TP53 R248W"
'
SNV_C1 <- SNV(
  chr = "17", 
  chr_pos = 7661779,  
  #alt = "C", 
  #ref = "G", 
  allele = 0, 
  #cause = ""
)
'

# STK11 LOH
CNA_C2 <- CNA(
  type = "D", 
  chr = "19", 
  chr_pos = 1177558, 
  len = 50873, 
  allele = 0, 
  #src_allele = NULL
)

# EGFR amp
CNA_C3 <- CNA(
  type = "A", 
  chr = "7", 
  chr_pos = 55019017,  
  len = 192611, 
  #allele = 2, 
  #src_allele = NULL
)

# KEAP1 R413C/H/L
SNV_C4 <- "KEAP1 R460M"
'
SNV_C4 <- SNV(
  chr = "19", 
  chr_pos = 10486125, 
  alt = "C", 
  ref = "G", 
  allele = 0, 
  #cause = "?"
)
'

# KRAS G12D
SNV_C5 <- "KRAS G12D"
'
SNV_C5 <- SNV(
  chr = "12", 
  chr_pos = 25205246, 
  #alt = "C", 
  #ref = "G", 
  allele = 0, 
  #cause = "?"
)
'

# KRAS amp (mutant)
CNA_C6 <- CNA(
  type = "A", 
  chr = "12", 
  chr_pos = 25205246, 
  len = 45690, 
  #allele = 2, 
  #src_allele = NULL
)


#CNA_C7 <- WGD


#-------------------------------------------------------------------------------
#-------------------------- add mutations to mutants ---------------------------
#-------------------------------------------------------------------------------

message("start placing mutations")
#m_engine$add_mutant(
#  mutant_name = "Clone 1", 
#  passenger_rates = c(SNV = mu_SNV), 
#  drivers = list(list("KRAS G12D", allele = 1))
#)

# Clone 1 : TP53 R248Q/W/L
m_engine$add_mutant(
  mutant_name = "C1", 
  passenger_rates = passengers, 
  driver_SNVs = list(SNV_C1)
)
message("C1 mutations added!")

# Clone 2 : STK11 LOH
m_engine$add_mutant(
  mutant_name = "C2", 
  passenger_rates = passengers, 
  driver_SNVs = list(CNA_C2)
)
message("C2 mutations added!")

# Clone 3 : EGFR amp
m_engine$add_mutant(
  mutant_name = "C3", 
  passenger_rates = passengers, 
  driver_SNVs = list(CNA_C3)
)
message("C3 mutations added!")

# Clone 4 : KEAP1 R413C/H/L
m_engine$add_mutant(
  mutant_name = "C4", 
  passenger_rates = passengers, 
  driver_SNVs = list(SNV_C4)
)
message("C4 mutations added!")

# Clone 5 : KRAS G12D
m_engine$add_mutant(
  mutant_name = "C5", 
  passenger_rates = passengers, 
  driver_SNVs = list(SNV_C5)
)
message("C5 mutations added")

# Clone 6 : KRAS amp (mutant)
m_engine$add_mutant(
  mutant_name = "C6", 
  passenger_rates = passengers, 
  driver_SNVs = list(CNA_C6)
)
message("C6 mutations added")

# Clone 6 : KRAS amp (mutant)
m_engine$add_mutant(
  mutant_name = "C7", 
  passenger_rates = passengers, 
  driver_SNVs = list(WGD)
)
message("C7 mutations added!")


#-------------------------------------------------------------------------------
#----------------------------- mutational signatures ---------------------------
#-------------------------------------------------------------------------------

timing <- readRDS( paste(curr_dir, "/data/chemo_timing.rds", sep = "") )

# SBS1, SBS4 and SBS5 are always active
m_engine$add_exposure(
  time = 0,
  coefficients = c(SBS1 = 0.35, SBS4 = 0.45, SBS5 = 0.2, ID1 = 1)
  )

m_engine$add_exposure(
  time = timing$chemo1_start, 
  coefficients = c(SBS1 = 0.35, SBS4 = 0.05, SBS5 = 0.15, SBS11 = 0.2, SBS25 = 0.25, ID1 = 1)
  )

m_engine$add_exposure(
  time = timing$chemo1_end, 
  coefficients = c(SBS1 = 0.4, SBS4 = 0.3, SBS5 = 0.1, SBS11 = 0.2, ID1 = 1)
)

m_engine$add_exposure(
  time = timing$chemo2_start, 
  coefficients = c(SBS1 = 0.4, SBS4 = 0.3, SBS5 = 0.1, SBS11 = 0.2, ID1 = 1)
)

m_engine$add_exposure(
  time = timing$chemo2_end, 
  coefficients = c(SBS1 = 0.4, SBS4 = 0.3, SBS5 = 0.1, SBS11 = 0.2, ID1 = 1)
)

#m_engine

#-------------------------------------------------------------------------------
#------------------------------ place mutations --------------------------------
#-------------------------------------------------------------------------------
samples_forest <- load_samples_forest( paste(curr_dir, "/data/forest.sff", sep = "") )

phylo_forest <- m_engine$place_mutations(samples_forest, num_of_preneoplatic_SNVs = 800, num_of_preneoplatic_indels = 200)

#-------------------------------------------------------------------------------
#------------------------------------- save ------------------------------------
#-------------------------------------------------------------------------------
phylo_forest$save( paste(curr_dir, "/data/phylo_forest.sff", sep = "") )

message("all mutations are placed!")

#-------------------------------------------------------------------------------
#----------------------------------- PLOTS -------------------------------------
#-------------------------------------------------------------------------------
annot_forest <- plot_forest(samples_forest) %>%
  annotate_forest(phylo_forest,
                  samples = T,
                  MRCAs = T,
                  exposures = T,
                  drivers=T,
                  add_driver_label = T)


p_exposure_timeline <- plot_exposure_timeline(phylo_forest)

labels <- get_relevant_branches(samples_forest)
sticks <- plot_sticks(samples_forest, labels)

pl <- annot_forest + sticks + p_exposure_timeline + plot_layout(nrow = 3, design = 'A\nA\nB\nB\nC')
#pl <- annot_forest + sticks + plot_layout(nrow = 3, design = 'A\nA\nB\nB\nC')
ggsave( paste(curr_dir, "/plots/SPN06_mutations.png", sep = "") , plot = pl, width = 210, height = 297, units = "mm", dpi = 300)


message("DONE!!!!!!!")




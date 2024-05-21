
library(rRACES)
library(dplyr)
library(ggplot2)

#-------------------------------------------------------------------------------
#------------------------ Customizing MutationEngine ---------------------------
#-------------------------------------------------------------------------------

# PRE-SET ENGINE
#-------------------------------------------------------------------------------
#m_engine <- build_mutation_engine(
#  setup_code = "GRCh38", context_sampling = 20
#)

# CUSTOMIZED ENGINE
#-------------------------------------------------------------------------------
reference_url <- paste0("https://ftp.ensembl.org/pub/grch37/release-111/",
                        "fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.",
                        "dna.chromosome.22.fa.gz")

sbs_url <- paste0("https://cancer.sanger.ac.uk/signatures/documents/2123/",
                  "COSMIC_v3.4_SBS_GRCh37.txt")

drivers_url <- paste0("https://raw.githubusercontent.com/",
                      "caravagnalab/rRACES/main/inst/extdata/",
                      "driver_mutations_hg19.csv")

passenger_cnas_url <- paste0("https://raw.githubusercontent.com/",
                             "caravagnalab/rRACES/main/inst/extdata/",
                             "passenger_CNAs_hg19.csv")

germline_url <- paste0("https://www.dropbox.com/scl/fi/g9oloxkip18tr1r",
                       "m6wjve/germline_data_demo.tar.gz?rlkey=15jshul",
                       "d3bqgyfcs7fa0bzqeo&dl=1")

m_engine <- build_mutation_engine(
  directory = "/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/engine", 
  reference_src = reference_url, 
  SBS_src = sbs_url, 
  drivers_src = drivers_url, 
  passenger_CNAs_src = passenger_cnas_url, 
  germline_src = germline_url
  )
dir.exists("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/engine")


plot_forest <- F

get_mutation_engine_codes()


#-------------------------------------------------------------------------------
#------------------------------ point mutations --------------------------------
#-------------------------------------------------------------------------------

# driver mutations
SNV_C1 <- SNV(chr = "22", chr_pos = 10510210, alt = "C", ref = "G", allele = 1, cause = "?")
SNV_C2 <- SNV(chr = "22", chr_pos = 10010210, alt = "C", ref = "G", allele = 1, cause = "?")
SNV_C3 <- SNV(chr = "22", chr_pos = 1051210,  alt = "C", ref = "G", allele = 1, cause = "?")
SNV_C4 <- SNV(chr = "22", chr_pos = 20510210, alt = "C", ref = "G", allele = 1, cause = "?")
SNV_C5 <- SNV(chr = "22", chr_pos = 100210,   alt = "C", ref = "G", allele = 1, cause = "?")
SNV_C6 <- SNV(chr = "22", chr_pos = 12510210, alt = "C", ref = "G", allele = 1, cause = "?")

passengers <- c(SNV = 2e-8, CNA = 1e-11)

# Clone 1 : TP53 R248Q/W/L
m_engine$add_mutant(
  mutant_name = "C1", 
  passenger_rates = passengers, 
  drivers = list(SNV_C1)
)

# Clone 2 : STK11 LOH
m_engine$add_mutant(
  mutant_name = "C2", 
  passenger_rates = passengers, 
  driver_SNVs = list(SNV_C2)
)

# Clone 3 : EGFR amp
m_engine$add_mutant(
  mutant_name = "C3", 
  passenger_rates = passengers, 
  driver_SNVs = list(SNV_C3)
)

# Clone 4 : KEAP1 R413C/H/L
m_engine$add_mutant(
  mutant_name = "C4", 
  passenger_rates = passengers, 
  driver_SNVs = list(SNV_C4)
)

# Clone 5 : KRAS G12D
m_engine$add_mutant(
  mutant_name = "C5", 
  passenger_rates = passengers, 
  driver_SNVs = list(SNV_C5)
)

# Clone 6 : KRAS amp (mutant)
m_engine$add_mutant(
  mutant_name = "C6", 
  passenger_rates = passengers, 
  driver_SNVs = list(SNV_C6)
)

#-------------------------------------------------------------------------------
#----------------------------- Mutational Signatures ---------------------------
#-------------------------------------------------------------------------------

timing <- readRDS("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/chemo_timing.rds")

# SBS1, SBS4 and SBS5 are always active
m_engine$add_exposure(
  coefficients = c(SBS1 = 0.5, SBS4 = 0.3, SBS5 = 0.2)
  )

m_engine$add_exposure(
  time = timing$chemo1_start, 
  coefficients = c(SBS1 = 0.35, SBS4 = 0.05, SBS5 = 0.15, SBS11 = 0.2, SBS25 = 0.25)
  )

m_engine$add_exposure(
  time = timing$chemo1_end, 
  coefficients = c(SBS1 = 0.4, SBS4 = 0.3, SBS5 = 0.1, SBS11 = 0.2)
)

m_engine

#-------------------------------------------------------------------------------
#------------------------------ Place Mutations --------------------------------
#-------------------------------------------------------------------------------

samples_forest <- load_samples_forest("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/forest.sff")

phylo_forest <- m_engine$place_mutations(samples_forest, 1000)



#-------------------------------------------------------------------------------
#------------------------------ SO FAR SO GOOD ---------------------------------
#-------------------------------------------------------------------------------

phylo_forest$get_exposures()

m_engine$get_SBSs()[1:6, 1:5]

plot_exposure_timeline(phylo_forest, linewidth = 0.8, emphatize_switches = FALSE)






all_SNV <- phylo_forest$get_sampled_cell_mutations() %>% as_tibble()
all_SNV %>%
  group_by(cause) %>%
  summarise(nPos = n_distinct(chr_pos))

all_SNV %>%
  group_by(class) %>%
  summarise(nPos = n_distinct(chr_pos))

all_CNV <- phylo_forest$get_sampled_cell_CNAs() %>% as_tibble()
all_CNV %>%
  group_by(class) %>%
  summarise(nPos = n_distinct(begin))

phylo_forest$save("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/phylo_forest.sff")

if (plot_forest == T) {
  annot_forest <- plot_forest(forest) %>%
    annotate_forest(phylo_forest, exposures = F)
  ggsave('/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/SPN06_ann_forest.png', plot = annot_forest,   height = 12, width = 10, dpi = 300, units = 'in')
}



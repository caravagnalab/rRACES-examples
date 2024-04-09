rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils.R")

seed <- 12345
set.seed(seed)

# Mutation generation ####
# mutation engine
# work on chr2
reference_url <- paste0("https://ftp.ensembl.org/pub/grch37/current/",
                        "fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.",
                        "dna.chromosome.2.fa.gz")

sbs_url <- paste0("https://cancer.sanger.ac.uk/signatures/documents/2123/",
                  "COSMIC_v3.4_SBS_GRCh37.txt")

drivers_url <- paste0("https://raw.githubusercontent.com/",
                      "caravagnalab/rRACES/main/inst/extdata/",
                      "driver_mutations_hg19.csv")

passenger_cnas_url <- paste0("https://raw.githubusercontent.com/",
                             "caravagnalab/rRACES/main/inst/extdata/",
                             "passenger_CNAs_hg19.csv")

germline_url <- paste0("https://www.dropbox.com/scl/fi/3rs2put4wde3objxmhvjc/germline_data_hg19.tar.gz?rlkey=imawitklf8d6zphz9ugriv4qm&dl=1")


# build a mutation engine and place all the files in the directory "Test"
get_mutation_engine_codes()
m_engine <- build_mutation_engine(setup_code = "GRCh37", context_sampling = 20)

m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = 5e-8, CNA = 1e-11),
                    driver_SNVs = c(SNV("2", 209113113, "A")))


m_engine$add_mutant("Clone 2",
                    passenger_rates = c(SNV=5e-8, CNA=1e-11),
                    driver_SNVs = c(),
                    driver_CNAs = c(CNA(type = "A", "6", pos_in_chr = 25100000,len = 1e7)))

m_engine$add_mutant("Clone 3",
                    passenger_rates = c(SNV=5e-8, CNA=1e-11),
                    driver_SNVs = c(SNV("1", 115256530, "T")),
                    driver_CNAs = c())

# Mutational exposures ####
treatment_info <- readRDS("data/treatment_info.rds")
m_engine$add_exposure(coefficients = c(SBS5 = 0.5, SBS1 = 0.5)) # ID1 is missing
m_engine$add_exposure(time = treatment_info$treatment_start, c(SBS5 = 0.4, SBS1 = 0.4, SBS25 = 0.2))
m_engine$add_exposure(time = treatment_info$treatment_end, coefficients = c(SBS5 = 0.5, SBS1 = 0.5))

# Sequence mutations ####
sampled_phylogeny <- load_samples_forest("data/samples_forest.sff")
phylo_forest <- m_engine$place_mutations(sampled_phylogeny, 1000)

tree_plot <- plot_forest(sampled_phylogeny)
tree_plot <- annotate_forest(tree_plot, forest = phylo_forest, exposures = TRUE, MRCAs = FALSE, samples = FALSE, facet_signatures = TRUE, drivers = TRUE, add_driver_label = FALSE)

phylo_forest$save("data/phylo_forest.sff")
ggsave("tissue/tree.pdf", dpi=300, width = 12, height = 12, plot = tree_plot)

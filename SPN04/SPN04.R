setwd("~/GitHub/rRACES-examples/SPN04")

rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils.R")

seed <- 777
set.seed(seed)
dt <- .1

# Prep simulation ####
sim <- new(Simulation, seed = seed, save_snapshot = F)
sim$update_tissue(1e3, 1e3)
# Set the "border" growth model
sim$duplicate_internal_cells <- TRUE

# Set the death activation level to avoid drift
sim$death_activation_level <- 50

# First and Second mutant ####
sim$add_mutant(name = "Clone 1", growth_rates = .5, death_rates = 0.01)

sim$place_cell("Clone 1", 500, 500)
# Let the simulation evolve until "A+" consists of 1000 cells
#sim$run_up_to_size("Clone 1", 1) # to observe muller from beginning
sim <- run_up_to_size_by_steps(sim, "Clone 1", 1000, dt)

sim$add_mutant("Clone 2",growth_rates = 1, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$get_species()
sim <- sweep_population(sim, "Clone 1", "Clone 2", .5, first_reduction = 10, delta_time = dt)

# Third mutant ####
sim$add_mutant("Clone 3",growth_rates = 2, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")

sim <- sweep_population(sim, "Clone 2", "Clone 3", .75, first_reduction = 10, delta_time = dt)
sim <- run_up_to_size_by_steps(sim, "Clone 3", .75e5, dt)

plot_tissue(sim)
ggsave("tissue/tissue_00.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim)
ggsave("tissue/muller_00.pdf", dpi=300, width = 8, height = 8)

# First sampling ####
n_w <- n_h <- 50
ncells <- .99*n_w*n_h
bbox <- sim$search_sample(c("Clone 3" = ncells), n_w, n_h)
sim$sample_cells("A", bbox$lower_corner, bbox$upper_corner)

plot_tissue(sim, num_of_bins = 300)
ggsave("tissue/tissue_01.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim)
ggsave("tissue/muller_01.pdf", dpi=300, width = 8, height = 8)

# Treatment ####
treatment_start <- sim$get_clock()
sim$update_rates("Clone 1",rates = c(growth = 0, death=5))
sim$update_rates("Clone 2",rates = c(growth = 0, death=5))
sim$update_rates("Clone 3",rates = c(growth = 0, death=5))
sim <- run_down_to_size_by_steps(sim, "Clone 3", 1e3, delta_time = .1)
treatment_end <- sim$get_clock()

plot_tissue(sim, num_of_bins = 300)
ggsave("tissue/tissue_02.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim) + xlim(20, NA)
ggsave("tissue/muller_02.pdf", dpi=300, width = 8, height = 8)

# Relapse ####

sim$update_rates("Clone 3",rates = c(growth = 2, death=0.01))
sim <- run_up_to_size_by_steps(sim, "Clone 3", .75e5, delta_time = .1)

plot_tissue(sim, num_of_bins = 300)
ggsave("tissue/tissue_03.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim) + xlim(20, NA)
ggsave("tissue/muller_03.pdf", dpi=300, width = 8, height = 8)

# Second sampling ####
n_w <- n_h <- 50
ncells <- .99*n_w*n_h
bbox <- sim$search_sample(c("Clone 3" = ncells), n_w, n_h)
sim$sample_cells("B", bbox$lower_corner, bbox$upper_corner)

plot_tissue(sim, num_of_bins = 300)
ggsave("tissue/tissue_04.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim) + xlim(20, NA)
ggsave("tissue/muller_04.pdf", dpi=300, width = 8, height = 8)

sampled_phylogeny <- sim$get_samples_forest()

# Recap Plot ####
muller_plot <- plot_muller(sim)
piechart <- plot_state(sim)
timeseries_plot <- plot_timeseries(sim)
final_sampled_tissue <- plot_tissue(sim)
plot_phylogeny = plot_forest(sampled_phylogeny) %>%
  annotate_forest(sampled_phylogeny, 
                  samples = TRUE, 
                  MRCAs = TRUE, 
                  exposures = FALSE, 
                  facet_signatures = FALSE, 
                  drivers = FALSE, 
                  add_driver_label = FALSE)
layout <- "
ACBD
#EE#
"

summary_patchwork = patchwork::wrap_plots(
  piechart, muller_plot, final_sampled_tissue, timeseries_plot,plot_phylogeny,
  guides = 'auto', design = layout
)
summary_patchwork

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
m_engine <- build_mutation_engine(directory = "Test_SPN04",
                                  reference_src = reference_url,
                                  SBS_src = sbs_url,
                                  drivers_src = drivers_url,
                                  passenger_CNAs_src = passenger_cnas_url,
                                  germline_src = germline_url)

m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = 1e-7, CNA = 0),
                    driver_SNVs = c(SNV("2", 209113113, "T")))
                    

m_engine$add_mutant("Clone 2",
                    passenger_rates = c(SNV=1e-7, CNA=1e-9),
                    driver_SNVs = c(SNV("2", 209113113, "T"))
                    #driver_SNVs = c(SNV("2", 209113113, "A"))#,
                    #driver_SNVs = c(SNV("2", 25457163, "T"))
                    #driver_CNAs = c(CNA(type = "A", "2", pos_in_chr = 209113113,len = 100))
)

m_engine$add_mutant("Clone 3",
                    passenger_rates = c(SNV=1e-7, CNA=1e-9),
                    driver_SNVs = c(SNV("2", 209113113, "T"))
                    #driver_SNVs = c(SNV("2", 209113113, "A"), SNV("2", 209113112, "T"))#,
                    #driver_SNVs = c(SNV("2", 25457163, "T"))
                    #driver_CNAs = c(CNA(type = "A", "2", pos_in_chr = 209113113,len = 100))
)

# Mutational exposures ####
m_engine$add_exposure(coefficients = c(SBS5 = 0.5, SBS1 = 0.5)) # ID1 is missing
m_engine$add_exposure(time = treatment_start, c(SBS5 = 0.4, SBS1 = 0.4, SBS25 = 0.2))
m_engine$add_exposure(time = treatment_end, coefficients = c(SBS5 = 0.5, SBS1 = 0.5))

# Sequence mutations ####
phylo_forest <- m_engine$place_mutations(sampled_phylogeny, 0)

tree_plot <- plot_forest(sampled_phylogeny)
annotate_forest(tree_plot, forest = NULL, exposures = TRUE)

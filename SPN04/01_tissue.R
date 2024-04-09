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
sim$add_mutant(name = "Clone 1", growth_rates = .05, death_rates = .001)

sim$place_cell("Clone 1", 500, 500)
# Let the simulation evolve until "A+" consists of 1000 cells
#sim$run_up_to_size("Clone 1", 1) # to observe muller from beginning
sim <- run_up_to_size_by_steps(sim, "Clone 1", 1000, dt)

sim$add_mutant("Clone 2",growth_rates = .1, death_rates = .001)
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$get_species()
sim <- sweep_population(sim, "Clone 1", "Clone 2", .5, first_reduction = 10, delta_time = dt)

# Third mutant ####
sim$add_mutant("Clone 3",growth_rates = .2, death_rates = .001)
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")

sim <- sweep_population(sim, "Clone 2", "Clone 3", .75, first_reduction = 10, delta_time = dt)
sim <- run_up_to_size_by_steps(sim, "Clone 3", .75e5, dt)

plot_tissue(sim)
ggsave("tissue/tissue_00.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim)
ggsave("tissue/muller_00.pdf", dpi=300, width = 8, height = 8)

# First sampling ####
n_w <- n_h <- 15
ncells <- .99*n_w*n_h
bbox <- sim$search_sample(c("Clone 3" = ncells), n_w, n_h)
sim$sample_cells("A", bbox$lower_corner, bbox$upper_corner)

plot_tissue(sim, num_of_bins = 300)
ggsave("tissue/tissue_01.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim)
ggsave("tissue/muller_01.pdf", dpi=300, width = 8, height = 8)

# Treatment ####
treatment_start <- sim$get_clock()
sim$update_rates("Clone 1",rates = c(growth = 0, death=.5))
sim$update_rates("Clone 2",rates = c(growth = 0, death=.5))
sim$update_rates("Clone 3",rates = c(growth = 0, death=.5))
sim <- run_down_to_size_by_steps(sim, "Clone 3", 1e3, delta_time = .1)
treatment_end <- sim$get_clock()

plot_tissue(sim, num_of_bins = 300)
ggsave("tissue/tissue_02.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim) + xlim(20, NA)
ggsave("tissue/muller_02.pdf", dpi=300, width = 8, height = 8)

# Relapse ####

sim$update_rates("Clone 3",rates = c(growth = .2, death=.001))
sim <- run_up_to_size_by_steps(sim, "Clone 3", .75e5, delta_time = .1)

plot_tissue(sim, num_of_bins = 300)
ggsave("tissue/tissue_03.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim) + xlim(20, NA)
ggsave("tissue/muller_03.pdf", dpi=300, width = 8, height = 8)

# Second sampling ####
n_w <- n_h <- 15
ncells <- .99*n_w*n_h
bbox <- sim$search_sample(c("Clone 3" = ncells), n_w, n_h)
sim$sample_cells("B", bbox$lower_corner, bbox$upper_corner)

plot_tissue(sim, num_of_bins = 300)
ggsave("tissue/tissue_04.pdf", dpi=300, width = 8, height = 8)
plot_muller(sim) + xlim(20, NA)
ggsave("tissue/muller_04.pdf", dpi=300, width = 8, height = 8)

sampled_phylogeny <- sim$get_samples_forest()
sampled_phylogeny$save("data/samples_forest.sff")
treatment_info <- list(treatment_start=treatment_start, treatment_end=treatment_end)
saveRDS(treatment_info, "data/treatment_info.rds")

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



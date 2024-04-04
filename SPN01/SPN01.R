rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
# Set directories
dir.create(path = "tissue", recursive = TRUE)
# Prep simulation ####
sim <- new(Simulation, seed = 5, save_snapshot = F)
sim$update_tissue(1e3, 1e3)
# Set the "border" growth model
sim$duplicate_internal_cells <- FALSE

# Set the death activation level to avoid drift
sim$death_activation_level <- 50

# Add first mutant
sim$add_mutant(name = "Clone 1",
               growth_rates = 1,
               death_rates = 0.01)

sim$place_cell("Clone 1", 500, 500)

# Let the simulation evolve until "Clone 1" consists of 1300 cells
sim$run_up_to_size("Clone 1", 1300)
clone1_tissue <-plot_tissue(sim)
clone1_muller <- plot_muller(sim)

# Add second mutant
sim$add_mutant("Clone 2",growth_rates = 2,
               death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$update_rates("Clone 1",rates = c(growth = 0, death=0))


# Let the simulation evolve until "Clone 2" consists of 5000 cells
sim$run_up_to_size("Clone 2", 100)
sim$update_rates("Clone 1",rates = c(growth = 1, death=0))
sim$run_up_to_size("Clone 2",5000)


clone2_tissue <-plot_tissue(sim)
clone2_muller <- plot_muller(sim)
ggsave("tissue/muller_02.pdf", clone2_muller,dpi=300, width = 8, height = 8)

# Add third mutant
sim$add_mutant("Clone 3",growth_rates = 5,
               death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$update_rates("Clone 2",rates = c(growth = 0, death=0))
sim$run_up_to_size("Clone 3", 100)
sim$update_rates("Clone 2",rates = c(growth = 2, death=0))
sim$run_up_to_size("Clone 3", 1e4)

clone3_tissue <-plot_tissue(sim)
clone3_muller <- plot_muller(sim)
ggsave("tissue/muller_03.pdf", clone3_muller,dpi=300, width = 8, height = 8)
# Add last mutant
sim$add_mutant("Clone 4",growth_rates = 15,
               death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 3"), "Clone 4")
sim$update_rates("Clone 3",rates = c(growth = 0, death=0))
sim$run_up_to_size("Clone 4", 100)
sim$update_rates("Clone 3",rates = c(growth = 5, death=0.01))

sim$run_up_to_size("Clone 4", 3e5)

clone4_tissue <-plot_tissue(sim)
clone4_muller <- plot_muller(sim)
ggsave("tissue/muller_04.pdf", clone4_muller,dpi=300, width = 8, height = 8)


## Sampling
# Three boxes with 'ncells' cells each
n_w <- n_h <- 50
ncells <- 0.8*n_w*n_h
bbox <- sim$search_sample(c("Clone 4" = ncells), n_w, n_h)

plot_tissue(sim) +
  geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
            ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
            fill = NA, color = "black")
sim$sample_cells("A", bbox$lower_corner, bbox$upper_corner)
plot_tissue(sim)
bbox <- sim$search_sample(c("Clone 4" = ncells), n_w, n_h)
plot_tissue(sim) +
  geom_rect(xmin = bbox$lower_corner[1]+100, xmax = bbox$upper_corner[1]+100,
            ymin = bbox$lower_corner[2]+100, ymax = bbox$upper_corner[2]+100,
            fill = NA, color = "black")
sim$sample_cells("B", bbox$lower_corner+100, bbox$upper_corner+100)
plot_tissue(sim)
bbox <- sim$search_sample(c("Clone 4" = ncells), n_w, n_h)
plot_tissue(sim) +
  geom_rect(xmin = bbox$lower_corner[1]+200, xmax = bbox$upper_corner[1]+200,
            ymin = bbox$lower_corner[2]+200, ymax = bbox$upper_corner[2]+200,
            fill = NA, color = "black")
sim$sample_cells("C", bbox$lower_corner+200, bbox$upper_corner+200)
# Get forest ####
sampled_phylogeny <- sim$get_samples_forest()
# Recap Plot ####
final_sampled_tissue <-plot_tissue(sim)
muller_plot <- plot_muller(sim)
piechart <- plot_state(sim)
timeseries_plot <- plot_timeseries(sim)
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
ggsave("tissue/summary_plot.pdf",summary_patchwork)

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
# Set initial growth rates
Clone1_gr <- 1
Clone2_gr <- 2
Clone3_gr <- 5
Clone4_gr <- 15
# Add first mutant
sim$add_mutant(name = "Clone 1",
               growth_rates = Clone1_gr,
               death_rates = 0.01)

sim$place_cell("Clone 1", 500, 500)

# Let the simulation evolve until "Clone 1" consists of 1300 cells
sim$run_up_to_size("Clone 1", 1300)
clone1_tissue <-plot_tissue(sim)
clone1_muller <- plot_muller(sim)

# Add second mutant
sim$add_mutant("Clone 2",growth_rates = Clone2_gr,
               death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$update_rates("Clone 1",rates = c(growth = 0, death=0))


# Let the simulation evolve until "Clone 2" consists of 5000 cells
sim$run_up_to_size("Clone 2", 100)
# sim$update_rates("Clone 1",rates = c(growth = Clone1_gr, death=0))
sim$run_up_to_size("Clone 2",5000)


clone2_tissue <-plot_tissue(sim)
clone2_muller <- plot_muller(sim)
ggsave("tissue/muller_02.pdf", clone2_muller,dpi=300, width = 8, height = 8)
ggsave("tissue/tissue_02.pdf", clone2_tissue,dpi=300, width = 8, height = 8)
# Add third mutant
sim$add_mutant("Clone 3",growth_rates = Clone3_gr,
               death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$update_rates("Clone 2",rates = c(growth = 0, death=0))
sim$run_up_to_size("Clone 3", 100)
# sim$update_rates("Clone 2",rates = c(growth = Clone2_gr, death=0))
sim$run_up_to_size("Clone 3", 1e4)

clone3_tissue <-plot_tissue(sim)
clone3_muller <- plot_muller(sim)
ggsave("tissue/muller_03.pdf", clone3_muller,dpi=300, width = 8, height = 8)
ggsave("tissue/tissue_03.pdf", clone3_tissue,dpi=300, width = 8, height = 8)

# Add last mutant
sim$add_mutant("Clone 4",growth_rates = Clone4_gr,
               death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 3"), "Clone 4")
sim$update_rates("Clone 3",rates = c(growth = 0, death=0))
sim$run_up_to_size("Clone 4", 100)
# sim$update_rates("Clone 3",rates = c(growth = Clone3_gr, death=0.01))

sim$run_up_to_size("Clone 4", 4e5)

clone4_tissue <-plot_tissue(sim)
clone4_muller <- plot_muller(sim)
clone4_state <- plot_state(sim)
ggsave("tissue/muller_04.pdf", clone4_muller,dpi=300, width = 8, height = 8)
ggsave("tissue/tissue_04.pdf", clone4_tissue,dpi=300, width = 8, height = 8)


## Sampling
# Define three boxes for sampling in
# different position of the tissue
# 
bbox_1 <- new(TissueRectangle, c(145, 270), c(195, 320))
bbox_2 <- new(TissueRectangle, c(250, 600), c(300, 650))
bbox_3 <- new(TissueRectangle, c(600, 270), c(650, 320))
plot_tissue(sim) +
  geom_rect(xmin = bbox_1$lower_corner[1], xmax = bbox_1$upper_corner[1],
            ymin = bbox_1$lower_corner[2], ymax = bbox_1$upper_corner[2],
            fill = NA, color = "red")+
  geom_rect(xmin = bbox_2$lower_corner[1], xmax = bbox_2$upper_corner[1],
            ymin = bbox_2$lower_corner[2], ymax = bbox_2$upper_corner[2],
            fill = NA, color = "black") +
  geom_rect(xmin = bbox_3$lower_corner[1], xmax = bbox_3$upper_corner[1],
            ymin = bbox_3$lower_corner[2], ymax = bbox_3$upper_corner[2],
            fill = NA, color = "blue")
sim$sample_cells("A", bbox_1$lower_corner, bbox_1$upper_corner)
sim$sample_cells("B", bbox_2$lower_corner, bbox_2$upper_corner)
sim$sample_cells("C", bbox_3$lower_corner, bbox_3$upper_corner)
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

rm(list = ls())
library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)

set.seed(12345)
sim <- SpatialSimulation(name = 'SPN01', seed = 12345)
sim$history_delta <- 1 
sim$death_activation_level <- 50

# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = 0.08, death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1000)


# Clone 2
sim$add_mutant(name = "Clone 2", growth_rates = 0.3, death_rates = 0.01)
sim$update_rates(species = "Clone 1", rates = c(growth = 0.06, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size(species = 'Clone 2', 700)



sim$update_rates(species = "Clone 1", rates = c(growth = 0.02, death = 0.03))
sim$update_rates(species = "Clone 2", rates = c(growth = 0.4, death = 0.005))
sim$run_up_to_size(species = 'Clone 2', 6000)

# Grow 2
sim$update_rates(species = "Clone 1", rates = c(growth = 0.001, death = 0.3))
sim$run_up_to_size(species = 'Clone 2', 10000)
muller <- plot_muller(sim)
plot_state(sim)


# Clone 3
sim$add_mutant(name = "Clone 3", growth_rates = 2, death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.2, death = 0.02))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 8000)
muller <- plot_muller(sim)

# Clone 4
sim$add_mutant(name = "Clone 4", growth_rates = 5, death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.01, death = 0.02))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 4")
sim$run_up_to_size("Clone 4", 12000)
muller <- plot_muller(sim)

muller
plot_state(sim)


bbox <- bbox <- sim$search_sample(c("Clone 3" = 300), 25, 25)
sim$sample_cells("Sample_B", bbox$lower_corner, bbox$upper_corner)
t1.1<- plot_tissue(sim)


bbox <- bbox <- sim$search_sample(c("Clone 4" = 300), 25, 25)
sim$sample_cells("Sample_C", bbox$lower_corner, bbox$upper_corner)
t1.2<- plot_tissue(sim)

bbox <- bbox <- sim$search_sample(c("Clone 2" = 180, 'Clone 3' = 120), 25, 25)
sim$sample_cells("Sample_A", bbox$lower_corner, bbox$upper_corner)
t1.3<- plot_tissue(sim)

# Forest
forest <- sim$get_samples_forest()
forest$save("data/samples_forest_new_dynamics.sff")
plt_forest <- plot_forest(forest) %>%
  annotate_forest(forest)
plt_forest
piechart <- plot_state(sim)
timeseries <- plot_timeseries(sim)
# ggsave(filename = "~/Desktop/test_forest.png",plot = plt_forest, height = 12, width = 10, dpi = 300, units = 'in')
# # Final plot
pl <- t1.1 + t1.2 + t1.3 + piechart + timeseries + muller + plt_forest + plot_layout(design = 'ABC\nEFG\nHHH\nHHH\nHHH') 
pl
ggsave("plots/SPN01_tissue_new_dynamics.png", plot = pl, height = 12, width = 10, dpi = 300, units = 'in')

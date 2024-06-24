library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)

set.seed(12345)
sim <- SpatialSimulation(name = 'SPN03', seed = 12345)
sim$history_delta <- 1 
sim$death_activation_level <- 50

# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = 0.08, death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1000)

# Sample A 
bbox <- sim$search_sample(c("Clone 1" = 300), 25, 25)
sim$sample_cells("Sample A", bbox$lower_corner, bbox$upper_corner)
t1 <- plot_tissue(sim)

# Clone 2
sim$add_mutant(name = "Clone 2", growth_rates = 0.3, death_rates = 0.01)
sim$update_rates(species = "Clone 1", rates = c(growth = 0.06, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size(species = 'Clone 2', 700)

bbox <- sim$search_sample(c("Clone 1" = 100, 'Clone 2' = 50), 20, 20)
sim$sample_cells("Sample B", bbox$lower_corner, bbox$upper_corner)
t2 <- plot_tissue(sim) 

sim$update_rates(species = "Clone 1", rates = c(growth = 0.02, death = 0.03))
sim$update_rates(species = "Clone 2", rates = c(growth = 0.4, death = 0.005))
sim$run_up_to_size(species = 'Clone 2', 6000)

# Grow 2
sim$update_rates(species = "Clone 1", rates = c(growth = 0.001, death = 0.3))
sim$run_up_to_size(species = 'Clone 2', 16000)
muller <- plot_muller(sim)

# Sample C
bbox <- sim$search_sample(c('Clone 2' = 400), 22, 22)
sim$sample_cells("Sample C", bbox$lower_corner, bbox$upper_corner)
t3 <- plot_tissue(sim)

# Clone 3
sim$add_mutant(name = "Clone 3", growth_rates = 2, death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.2, death = 0.02))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 10000)
muller <- plot_muller(sim)

# Sample D
bbox <- sim$search_sample(c("Clone 2" = 200, 'Clone 3' = 150), 20, 20)
sim$sample_cells("Sample D", bbox$lower_corner, bbox$upper_corner)
t4 <- plot_tissue(sim)

# Forest
forest <- sim$get_samples_forest()
forest$save("/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/samples_forest.sff")

plt_forest <- plot_forest(forest) %>%
  annotate_forest(forest)

piechart <- plot_state(sim)
timeseries <- plot_timeseries(sim)

# Final plot
pl <- t1 + t2 + t3 + t4  + piechart + timeseries + muller + plt_forest + plot_layout(design = 'ABCD\nEFGG\nHHHH\nHHHH\nHHHH') 
ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/plots/SPN03_tissue.png', plot = pl, width = 210, height = 297, units = "mm", dpi = 300)

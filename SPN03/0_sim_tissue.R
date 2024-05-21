library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)

sim <- new(Simulation, seed = 12345)
sim$duplicate_internal_cells <- FALSE
sim$history_delta <- 1
sim$death_activation_level <- 50

# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = 0.1, death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1000)

# Sample A 
bbox <- sim$search_sample(c("Clone 1" = 200), 20, 20)
sim$sample_cells("Sample A", bbox$lower_corner, bbox$upper_corner)
t1 <- plot_tissue(sim)


# Clone 2
sim$add_mutant(name = "Clone 2", growth_rates = 0.2, death_rates = 0.01)
sim$update_rates(species = "Clone 1", rates = c(growth = 0.08, death = 0.05))
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size(species = 'Clone 2', 1500)


bbox <- sim$search_sample(c("Clone 1" = 180, 'Clone 2' = 20), 20, 20)
sim$sample_cells("Sample B", bbox$lower_corner, bbox$upper_corner)
t2 <- plot_tissue(sim) 

sim$update_rates(species = "Clone 1", rates = c(growth = 0.05, death = 0.08))
sim$run_up_to_size(species = 'Clone 2', 5000)


# Grow 2
sim$update_rates(species = "Clone 1", rates = c(growth = 0, death = 0.2))
sim$run_up_to_size(species = 'Clone 2', 16000)

# Sample C
bbox <- sim$search_sample(c('Clone 2' = 200), 20, 20)
sim$sample_cells("Sample C", bbox$lower_corner, bbox$upper_corner)
t3 <- plot_tissue(sim)

# Clone 3
sim$add_mutant(name = "Clone 3", growth_rates = 2, death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.2, death = 0.02))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 15000)
muller <- plot_muller(sim)
ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/new_FALSE_1_muller_SPN03.png', plot = muller, dpi = 300, units = 'in')



# Sample D
bbox <- sim$search_sample(c("Clone 2" = 120, 'Clone 3' = 100), 20, 20)
sim$sample_cells("Sample D", bbox$lower_corner, bbox$upper_corner)
t4 <- plot_tissue(sim)

# Forest
forest <- sim$get_samples_forest()
#forest$save("/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/samples_forest.sff")

plt_forest <- plot_forest(forest) %>%
  annotate_forest(forest)
ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/new_FALSE_1_SPN03_forest.png', plt_forest)


piechart <- plot_state(sim)
timeseries <- plot_timeseries(sim)

# Final plot
pl <- t1 + t2 + t3 + t4  + piechart + timeseries + muller + plt_forest + plot_layout(design = 'ABCD
                                                                                         EFGG
                                                                                         HHHH
                                                                                         HHHH
                                                                                         HHHH
                                                                                         HHHH')
ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/new_FALSE_1_SPN03.png', plot = pl, height = 12, width = 10, dpi = 300, units = 'in')

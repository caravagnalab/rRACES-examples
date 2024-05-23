set.seed(123)
library(rRACES)
library(dplyr)
library(ggplot2)

sim <- new(Simulation, "SPN07",
           seed = 3,
           save_snapshot = F)
sim$update_tissue("Brain", 2e3, 2e3)
sim$duplicate_internal_cells <- T
sim$history_delta <- 0.1

sim$add_mutant(name = "1",
               growth_rates = 1.3,
               death_rates = 0)
sim$place_cell("1", 1000, 1000)
sim$run_up_to_size("1",1e4) #1e2
p1 = plot_tissue(sim, num_of_bins=100)
ggsave(p1, filename='./plots/tissue1.png')

sim$add_mutant(name = "2",
               growth_rates = 4,
               death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("1"), "2")
sim$run_up_to_size("2",1e4) #1e2
p2 = plot_tissue(sim, num_of_bins=100)
ggsave(p2, filename='./plots/tissue2.png')


sim$add_mutant(name = "3",
               growth_rates = 6,
               death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("2"), "3")
sim$run_up_to_size("3",1e4) #1e2
p3 = plot_tissue(sim, num_of_bins=100)
ggsave(p3, filename='./plots/tissue3.png')


sim$add_mutant(name = "4",
               growth_rates = 8,
               death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("2"), "4")
sim$update_rates("1", c(death = 10))
sim$run_up_to_size("4",1.1e5)
p4 = plot_tissue(sim, num_of_bins=100)
ggsave(p4, filename='./plots/tissue4.png')

#tissue_pre = plot_tissue(sim,num_of_bins  = 100)
#tissue_pre

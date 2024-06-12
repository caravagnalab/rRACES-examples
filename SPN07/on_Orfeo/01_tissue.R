set.seed(58)
library(rRACES)
library(dplyr)
library(ggplot2)
#setwd("~/dati_Orfeo/SPN07/on_Orfeo")
source('utils.R')

sim <- new(Simulation, "SPN07",
           seed = 3,
           save_snapshot = F)
sim$update_tissue("Brain", 2e3, 2e3)
sim$duplicate_internal_cells <- T
sim$history_delta <- 0.1

sim$add_mutant(name = "1",
               growth_rates = .2,
               death_rates = .03)
sim$place_cell("1", 1000, 1000)
sim$run_up_to_size("1",1e2) #1e2
p1 = plot_tissue(sim)
ggsave(p1, filename='./plots/tissue1.png')

sim$add_mutant(name = "2",
               growth_rates = .4,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("1"), "2")
sim$run_up_to_size("2",1e2) #1e2
p2 = plot_tissue(sim, num_of_bins=100)
ggsave(p2, filename='./plots/tissue2.png')

sim$update_rates("1", c(death = 1))
sim$add_mutant(name = "3",
               growth_rates = .7,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("2"), "3")
sim$run_up_to_size("3",1e2) #1e2
p3 = plot_tissue(sim, num_of_bins=100)
ggsave(p3, filename='./plots/tissue3.png')

#my_muller_plot(sim)

sim$add_mutant(name = "4",
               growth_rates = .9,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("2"), "4")
#sim$update_rates("1", c(death = 1))
sim$run_up_to_size("4",1e5)
p4 = plot_tissue(sim, num_of_bins=100)
ggsave(p4, filename='./plots/tissue4.png')

m1 = my_muller_plot(sim)
ggsave(m1, filename='./plots/log_muller.png')
m2 = rRACES::plot_muller(sim)
ggsave(m2, filename='./plots/muller_plot1.png')

#tissue_pre = plot_tissue(sim,num_of_bins  = 100)
#tissue_pre

### Sampling 1
sample_a = sim$search_sample(c('3'= 500), 50,50)
sample_b = sim$search_sample(c('4'=300), 50, 50)
sample_c = sim$search_sample(c('3'=300,'2'=200), 50, 50)

sim$sample_cells("A", sample_a$lower_corner, sample_a$upper_corner)
sim$sample_cells("B", sample_b$lower_corner, sample_b$upper_corner)
sim$sample_cells("C", sample_c$lower_corner, sample_c$upper_corner)

forest <- sim$get_samples_forest()
f1 = annotate_forest(forest = forest,tree_plot = plot_forest(forest),MRCAs = T)
ggsave(f1, filename='./plots/forest1.png')


### Chemotherapy
sim$update_rates("2", c(death = 5))
sim$update_rates("3", c(death = 5))
sim$update_rates("4", c(death = 5))
sim$run_up_to_time(sim$get_clock() + 2)
p5 = plot_tissue(sim)
ggsave(p5, filename='./plots/after_chemotherapy.png')

### Resistant clones birth
sim$add_mutant(name = "5",
               growth_rates = .3,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("4"), "5")
sim$run_up_to_size("5",1e3)

sim$add_mutant(name = "6",
               growth_rates = .5,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("5"), "6")
#sim$update_rates("4", c(death = 1e5))
sim$run_up_to_size("6",1e4)

p6 = plot_tissue(sim)
ggsave(p6, filename='./plots/tissue4.png')

m2 = my_muller_plot(sim)
ggsave(m2, filename='./plots/log_muller2.png')
m2 = rRACES::plot_muller(sim)
ggsave(m2, filename='./plots/muller_plot2.png')

### Sampling 2
sample_d = sim$search_sample(c('5'= 100), 50,50)
sample_e = sim$search_sample(c('6'=100), 50, 50)

sim$sample_cells("D", sample_d$lower_corner, sample_d$upper_corner)
sim$sample_cells("E", sample_e$lower_corner, sample_e$upper_corner)

forest <- sim$get_samples_forest()
f2 = annotate_forest(forest = forest,tree_plot = plot_forest(forest),MRCAs = T)
ggsave(f2, filename='./plots/forest2.png')
forest$save('forest_sampling_2.sff')



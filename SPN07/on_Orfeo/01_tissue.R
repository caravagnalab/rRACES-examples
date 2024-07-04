set.seed(58)
library(rRACES)
library(dplyr)
library(ggplot2)
#setwd("~/dati_Orfeo/orfeo/cephfs/home/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo")
#source('utils.R')

sim <- SpatialSimulation("SPN07",
           #seed = 3,
           save_snapshot = F,
           width = 2e3, height = 2e3)
#sim$update_tissue("Brain")
#sim$duplicate_internal_cells <- T
sim$border_growth_model = F
sim$history_delta <- 0.1

sim$add_mutant(name = "1",
               growth_rates = .1,
               death_rates = .03)
sim$place_cell("1", 1000, 1000)
sim$run_up_to_size("1",10) #1e2
p1 = plot_tissue(sim)
ggsave(p1, filename='./plots/tissue1.png')

sim$add_mutant(name = "2",
               growth_rates = .2,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("1"), "2")
sim$run_up_to_size("2",1e2) #1e2
p2 = plot_tissue(sim, num_of_bins=100)
ggsave(p2, filename='./plots/tissue2.png')

sim$update_rates("1", c(death = 1))
sim$add_mutant(name = "3",
               growth_rates = .3,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("2"), "3")
sim$run_up_to_size("3",1e2) #1e2
p3 = plot_tissue(sim, num_of_bins=100)
ggsave(p3, filename='./plots/tissue3.png')

#my_muller_plot(sim)

sim$add_mutant(name = "4",
               growth_rates = .5,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("2"), "4")
#sim$update_rates("1", c(death = 1))
sim$run_up_to_size("4",1e5)
p4 = plot_tissue(sim, num_of_bins=100)
ggsave(p4, filename='./plots/tissue4.png')

#m1 = my_muller_plot(sim)
#ggsave(m1, filename='./plots/log_muller.png')
m2 = rRACES::plot_muller(sim)
ggsave(m2, filename='./plots/muller_plot1.png')

state1= plot_state(sim)

#tissue_pre = plot_tissue(sim,num_of_bins  = 100)
#tissue_pre

### Sampling 1
sample_a = sim$search_sample(c('3'= 100), 25,25)
sample_b = sim$search_sample(c('4'=100), 25, 25)
sample_c = sim$search_sample(c('3'=50,'2'=50), 25, 25)

# plot_tissue(sim, num_of_bins=100) +
#   geom_rect(aes(xmin=sample_a$lower_corner[1],
#                 xmax=sample_a$upper_corner[1],
#                 ymin=sample_a$lower_corner[2],
#                 ymax=sample_a$upper_corner[2]),
#             color= 'white', fill='white')+
#   geom_rect(aes(xmin=sample_b$lower_corner[1],
#                 xmax=sample_b$upper_corner[1],
#                 ymin=sample_b$lower_corner[2],
#                 ymax=sample_b$upper_corner[2]),
#             color= 'white', fill='white')+
#   geom_rect(aes(xmin=sample_c$lower_corner[1],
#                 xmax=sample_c$upper_corner[1],
#                 ymin=sample_c$lower_corner[2],
#                 ymax=sample_c$upper_corner[2]),
#             color= 'white', fill='white')

sim$sample_cells("A", sample_a$lower_corner, sample_a$upper_corner)
sim$sample_cells("B", sample_b$lower_corner, sample_b$upper_corner)
sim$sample_cells("C", sample_c$lower_corner, sample_c$upper_corner)

forest <- sim$get_samples_forest()
f1 = annotate_forest(forest = forest,tree_plot = plot_forest(forest),MRCAs = T)
#ggsave(f1, filename='./plots/forest1.png')


### Chemotherapy
sim$update_rates("2", c(death = 3))
sim$update_rates("3", c(death = 3))
sim$update_rates("4", c(death = 2))
sim$run_up_to_time(sim$get_clock() + 2)
p5 = plot_tissue(sim)
ggsave(p5, filename='./plots/after_chemotherapy.png')

### Resistant clones birth
sim$add_mutant(name = "5",
               growth_rates = .3,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("4"), "5")
sim$run_up_to_size("5",1e2)

sim$add_mutant(name = "6",
               growth_rates = .5,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("5"), "6")
#sim$update_rates("4", c(death = 1e5))
sim$run_up_to_size("6",1e5)

p6 = plot_tissue(sim)
ggsave(p6, filename='./plots/tissue4.png')
state2= plot_state(sim)

#m2 = my_muller_plot(sim)
#ggsave(m2, filename='./plots/log_muller2.png')
m2 = rRACES::plot_muller(sim)
ggsave(m2, filename='./plots/muller_plot2.png')

### Sampling 2
sample_d = sim$search_sample(c('5'= 100), 25,25)
sample_e = sim$search_sample(c('6'=100), 25, 25)

sim$sample_cells("D", sample_d$lower_corner, sample_d$upper_corner)
sim$sample_cells("E", sample_e$lower_corner, sample_e$upper_corner)

forest <- sim$get_samples_forest()
f2 = annotate_forest(forest = forest,tree_plot = plot_forest(forest),MRCAs = T)
ggsave(f2, filename='./plots/forest2.png')
forest$save('forest_sampling.sff')
time_series = plot_timeseries(sim)


st = 'ABCDEF
      GHIILL
      MMMOOO
      MMMOOO'
report = patchwork::wrap_plots(p1,p2,p3,p4,p5,p6,
                               state1, state2, time_series,m2,
                               f1, f2,
                               design = st)

ggsave(report, filename= './plots/tissue_report.pdf',width = 20, height = 20)








# devtools::install_github("caravagnalab/rRACES")


library(rRACES)
library(dplyr)

sim <- new(Simulation, "SPN07",
           seed = 3,
           save_snapshot = F)

sim$duplicate_internal_cells <- T

sim$update_tissue("Liver", 2e3, 2e3)


sim$add_mutant(name = "1",
               growth_rates = 2,
               death_rates = 0)


sim$place_cell("1", 1000, 1000)

sim$run_up_to_size("1",1e2)

sim$add_mutant(name = "2",
               growth_rates = 4,
               death_rates = 0)

sim$mutate_progeny(sim$choose_cell_in("1"), "2")

sim$run_up_to_size("2",1e2)

print(sim)

sim$add_mutant(name = "3",
               growth_rates = 6,
               death_rates = 0)


sim$mutate_progeny(sim$choose_cell_in("2"), "3")


sim$run_up_to_size("3",1e2)

print(sim)

sim$add_mutant(name = "4",
               growth_rates = 8,
               death_rates = 0)

sim$mutate_progeny(sim$choose_cell_in("2"), "4")


sim$run_up_to_size("4",1e5)

print(sim)


tissue_pre = plot_tissue(sim,num_of_bins  = 50)

# sim = recover_simulation("SPN07")


bbox = sim$search_sample(c("3" = 100),nw = 50,nh = 50)


# bbox = tibble(lower_corner = c(750,1000), upper_corner = c(800,1080))


# plot_tissue(sim,num_of_bins  = 50)  +
#   geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
#             ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
#             fill = NA, color = "black")


sim$sample_cells("A", bbox$lower_corner, bbox$upper_corner)


forest <- sim$get_samples_forest()

# forest$get_nodes() %>% filter(!is.na(sample)) %>% group_by(mutant) %>% summarize(n = length(cell_id))
# 
# plot_forest(forest)


# bbox = tibble(lower_corner = c(1300,1060), upper_corner = c(1400,1080))

bbox = sim$search_sample(c("2" = 100,"3" = 100),nw = 50,nh = 50)


# plot_tissue(sim,num_of_bins  = 50)  +
#   geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
#             ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
#             fill = NA, color = "black")



sim$sample_cells("B", bbox$lower_corner, bbox$upper_corner)


forest <- sim$get_samples_forest()

# annotate_forest(tree_plot = plot_forest(forest),forest = forest)
# 
# forest$get_nodes() %>% filter(!is.na(sample),sample == "B") %>% 
#   group_by(mutant) %>% summarize(n = length(cell_id))


# bbox = sim$search_sample(c("3" = 300),nw = 50,nh = 50)

bbox = sim$search_sample(c("2" = 100,"4" = 100),nw = 50,nh = 50)


sim$sample_cells("C", bbox$lower_corner, bbox$upper_corner)


forest <- sim$get_samples_forest()

# annotate_forest(tree_plot = plot_forest(forest),forest = forest)
# 
# forest$get_nodes() %>% filter(!is.na(sample),sample == "C") %>% 
#   group_by(mutant) %>% summarize(n = length(cell_id))


sim$update_rates("1", c(death = 10))
sim$update_rates("2", c(death = 10))
sim$update_rates("3", c(death = 10))
sim$update_rates("4", c(death = 10))

sim$run_up_to_time(sim$get_clock() + 2)

print(sim)


sim$add_mutant(name = "5",
               growth_rates = 8,
               death_rates = 0)

sim$mutate_progeny(as.data.frame((sim$get_cells()  %>% 
              filter(mutant == '4') %>% arrange(position_x))[1,]), "5")

sim$run_up_to_size("5",1e2)

tissue_post = plot_tissue(sim)

sim$add_mutant(name = "6",
               growth_rates = 10,
               death_rates = 0)

sim$mutate_progeny(sim$choose_cell_in("5"), "6")

sim$run_up_to_size("6",1e5)

print(sim)

tissue_rel = plot_tissue(sim)

time_series = plot_timeseries(sim)

muller = plot_muller(sim)

# sim = recover_simulation("SPN07")

bbox = sim$search_sample(c("5" = 100),nw = 100,nh = 100)

# bbox = tibble(lower_corner = c(1200,1200), upper_corner = c(1260,1260))
# 
# 
# plot_tissue(sim,num_of_bins  = 50)  +
#   geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
#             ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
#             fill = NA, color = "black")

sim$sample_cells("D", bbox$lower_corner, bbox$upper_corner)

forest <- sim$get_samples_forest()

# annotate_forest(tree_plot = plot_forest(forest),forest = forest)

bbox = sim$search_sample(c("6" = 100),nw = 70,nh = 70)

sim$sample_cells("E", bbox$lower_corner, bbox$upper_corner)

forest <- sim$get_samples_forest()

plot_forest = annotate_forest(tree_plot = plot_forest(forest),forest = forest)


tissue_pre

tissue_post

tissue_rel

plot_forest

time_series 

muller



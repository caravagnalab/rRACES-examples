# devtools::install_github("caravagnalab/rRACES")
setwd("~/Documents/RACES")
library(rRACES)

sim <- new(Simulation, "SPN07",
           seed = 3,
           save_snapshot = F)

sim$duplicate_internal_cells <- T
sim$update_tissue("Brain", 2e3, 2e3)

# Clone with PTEN R
sim$add_mutant(name = "Clone_1", # A
               growth_rates = 2,
               death_rates = 0)

sim$place_cell("Clone_1", 1000, 1000)
sim$run_up_to_size("Clone_1",1e2)
# Clone with PTEN LOH
sim$add_mutant(name = "Clone_2", # NF1
               growth_rates = 10,
               death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("Clone_1"), "Clone_2")
sim$run_up_to_size("Clone_2",1e3)
print(sim)
plot_tissue(sim)
# Clone with NF1
sim$add_mutant(name = "Clone_3",
               growth_rates = 32,
               death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("Clone_2"), "Clone_3")
sim$run_up_to_size("Clone_3",1e4)
# print(sim)
plot_tissue(sim)
# Clone with ATRX
sim$add_mutant(name = "Clone_4",
               growth_rates = 120,
               death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("Clone_2"), "Clone_4")
sim$run_up_to_size("Clone_4",1e5)
#print(sim)
ggsave(plot_tissue(sim),filename='SamplingPoint1.png')

# Sample A
# bbox = sim$search_sample(c('Clone_2'=100,'Clone_3'=900) ,
#                          nw = 50,nh = 50
#                          )
# sim = recover_simulation("SPN07")
bbox = tibble(lower_corner = c(1000,1000), upper_corner = c(1300,1300))

plot_tissue(sim,num_of_bins  = 50)  +
  geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
            ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
            fill = NA, color = "black")
sim$sample_cells("Sample_A", bbox$lower_corner, bbox$upper_corner)

forest_sampleA <- sim$get_samples_forest()
plot_forest(forest_sampleA)

# forest$get_nodes() %>% filter(!is.na(sample)) %>% group_by(mutant) %>% summarize(n = length(cell_id))


bbox = tibble(lower_corner = c(1300,1060), upper_corner = c(1400,1080))


plot_tissue(sim,num_of_bins  = 50)  +
  geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
            ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
            fill = NA, color = "black")




sim$sample_cells("B", bbox$lower_corner, bbox$upper_corner)


forest <- sim$get_samples_forest()

plot_forest(forest)

forest$get_nodes() %>% filter(!is.na(sample),sample == "B") %>% 
  group_by(mutant) %>% summarize(n = length(cell_id))




sim$update_rates("A", c(death = 125))
sim$update_rates("B", c(death = 125))
sim$update_rates("C", c(death = 125))
sim$update_rates("D", c(death = 125))

sim$run_up_to_time(sim$get_clock() + 1)


sim$add_mutant(name = "E",
               growth_rates = 120,
               death_rates = 0)

sim$mutate_progeny(as.data.frame((sim$get_cells()  %>% 
                                    filter(mutant == 'D') %>% arrange(position_x))[1,]), "E")

sim$run_up_to_size("E",1e4)

plot_tissue(sim)

sim$add_mutant(name = "F",
               growth_rates = 400,
               death_rates = 0)

sim$mutate_progeny(sim$choose_cell_in("E"), "F")

sim$run_up_to_size("F",3e5)

print(sim)

plot_tissue(sim)

# sim = recover_simulation("SPN07")

bbox = sim$search_sample(c("E" = 1000),nw = 50,nh = 50)


plot_tissue(sim,num_of_bins  = 50)  +
  geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
            ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
            fill = NA, color = "black")

sim$sample_cells("C", bbox$lower_corner, bbox$upper_corner)


forest <- sim$get_samples_forest()

plot_forest(forest)


bbox = sim$search_sample(c("F" = 500,"E" = 500),nw = 50,nh = 50)

plot_tissue(sim,num_of_bins  = 50)  +
  geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
            ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
            fill = NA, color = "black")


sim$sample_cells("D", bbox$lower_corner, bbox$upper_corner)

forest <- sim$get_samples_forest()

plot_forest(forest)
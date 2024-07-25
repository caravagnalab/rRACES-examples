rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seed = 12345
set.seed(seed)

# Prep simulation ####
sim <- SpatialSimulation("SPN01", width=2e3, height=2e3)
#sim <- new(Simulation, seed = seed, save_snapshot = F)
sim$history_delta <- 1
#sim$update_tissue(2e3, 2e3)

# Set the death activation level to avoid drift
sim$death_activation_level <- 50

# First and Second mutant ####
sim$add_mutant(name = "Clone 1", growth_rates = .1, death_rates = .01)

####### CLONE 1 #######

sim$place_cell("Clone 1", 1000, 1000)
sim$run_up_to_size("Clone 1", 1000)
# plot_muller(sim)
t1 <- plot_tissue(sim,num_of_bins = 300)
####### CLONE 2 #######

sim$add_mutant("Clone 2",growth_rates = .3, death_rates = .01)
sim$update_rates("Clone 1", rates = c(growth = .033))
# plot_muller(sim)
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size("Clone 2", 1000)
# plot_muller(sim)

sim$update_rates("Clone 1", rates = c(growth = 0, death = 0.025))
sim$run_up_to_size("Clone 2", 2000)
plot_muller(sim)
# 
sim$update_rates("Clone 1", rates = c(growth = 0, death = .05))
sim$run_up_to_size("Clone 2", 4000)

t2 <- plot_tissue(sim,num_of_bins = 300)



####### CLONE 3 #######

sim$add_mutant("Clone 3",growth_rates = 1, death_rates = .01)
#sim$update_rates("Clone 2", rates = c(growth = .1, death = .05))
sim$update_rates("Clone 2", rates = c(growth = .1))

sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 4000)
plot_muller(sim)

sim$update_rates("Clone 2", rates = c(growth = 0, death = .025))
sim$run_up_to_size("Clone 3", 6000)
plot_muller(sim)

sim$update_rates("Clone 2", rates = c(growth = 0, death = .05))
sim$run_up_to_size("Clone 3", 1e5)
plot_muller(sim)

# sim$update_rates("Clone 2", rates = c(growth = 0, death = .075))
# sim$run_up_to_size("Clone 3", 1e5)
# plot_muller(sim)
t3 <- plot_tissue(sim,num_of_bins = 300)
####### CLONE 4 #######

sim$add_mutant("Clone 4",growth_rates = 2, death_rates = .01)
#sim$update_rates("Clone 2", rates = c(growth = .1, death = .05))
sim$update_rates("Clone 3", rates = c(growth = .1))
sim$mutate_progeny(sim$choose_cell_in("Clone 3"), "Clone 4")
sim$run_up_to_size("Clone 4", 1.5e5)
plot_muller(sim)

sim$update_rates("Clone 3", rates = c(growth = 0, death = .025))
sim$run_up_to_size("Clone 4", 2e5)
plot_muller(sim)

sim$update_rates("Clone 3", rates = c(growth = 0, death = 0.05))
sim$run_up_to_size("Clone 4", 2.5e5)
plot_muller(sim)

sim$update_rates("Clone 3", rates = c(growth = 0, death = 0.08))
sim$run_up_to_size("Clone 4", 3e5)
plot_muller(sim)

sim$update_rates("Clone 3", rates = c(growth = 0, death = 0.1))
sim$run_up_to_size("Clone 4", 3.5e5)
plot_muller(sim)

sim$update_rates("Clone 3", rates = c(growth = 0, death = 1))
sim$run_up_to_size("Clone 4", 3.8e5)
plot_muller(sim)

sim$update_rates("Clone 3", rates = c(growth = 0, death =100))
sim$update_rates("Clone 1", rates = c(growth = 0, death =100))
sim$update_rates("Clone 2", rates = c(growth = 0, death =100))

sim$run_up_to_size("Clone 4", 4e5)

muller_final<-plot_muller(sim)
t4<-plot_tissue(sim,num_of_bins = 300)
timeseries_final <- plot_timeseries(sim)
state_final <- plot_state(sim)


## Perform manual sampling to allow samples to be 
## as far as possible

t4_sampled <- t4 +
  ggplot2::geom_rect(xmin = 1100, xmax = 1119, ymin=600, ymax=619,fill = NA, color = "black") +
  ggplot2::geom_rect(xmin = 739, xmax = 758, ymin=1123, ymax=1142,fill = NA, color = "black") +
  ggplot2::geom_rect(xmin = 550, xmax = 569, ymin=550, ymax=569,fill = NA, color = "black")


sim$sample_cells("SPN01_Sample_1", lower_corner=c(1100,600), upper_corner=c(1119,619))
sim$sample_cells("SPN01_Sample_2", lower_corner=c(739,1123), upper_corner=c(758,1142))
sim$sample_cells("SPN01_Sample_3", lower_corner=c(550,550), upper_corner=c(569,569))



#n_w <- n_h <- 20
#ncells <- 0.8 * n_w * n_h
#
## Sampling ncells with random box sampling of boxes of size n_w x n_h
#bboxes <- sim$search_samples(c("Clone 4" = ncells), n_w, n_h,n_samples=3)
#counter <- 1
#for (bbox in bboxes) {
#  sample_name <-paste0("SPN01_Sample_",counter)
#  plot <- plot +
#    ggplot2::geom_rect(xmin = bbox$lower_corner[1],
#                       xmax = bbox$upper_corner[1],
#                       ymin = bbox$lower_corner[2],
#                       ymax = bbox$upper_corner[2],
#                       fill = NA, color = "black")
#  sim$sample_cells(sample_name, bbox$lower_corner, bbox$upper_corner)
#  counter <- counter+1
#}
#
forest <- sim$get_samples_forest()
forest$save("data/samples_forest.sff")
forest_final<- plot_forest(forest) %>%
  annotate_forest(forest,samples = T,MRCAs = T)


# Final plot
pl <- t1 +t2 +t3 +t4_sampled  +state_final + timeseries_final + muller_final + forest_final + plot_layout(design = 'ABCD
                                                                                         EFGG
                                                                                         HHHH
                                                                                         HHHH
                                                                                         HHHH
                                                                                         HHHH')
pl
ggsave("plots/SPN01_tissue.png", plot = pl, height = 12, width = 10, dpi = 300, units = 'in')

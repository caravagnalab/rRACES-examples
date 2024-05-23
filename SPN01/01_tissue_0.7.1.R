rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

seed = 12345
set.seed(seed)

# Prep simulation ####
sim <- new(Simulation, seed = seed, save_snapshot = F)
sim$history_delta <- 1
sim$update_tissue(2e3, 2e3)

# Set the death activation level to avoid drift
sim$death_activation_level <- 50

# First and Second mutant ####
sim$add_mutant(name = "Clone 1", growth_rates = .1, death_rates = .01)

####### CLONE 1 #######

sim$place_cell("Clone 1", 1000, 1000)
sim$run_up_to_size("Clone 1", 1000)
plot_muller(sim)
plot_tissue(sim)
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

plot_muller(sim)



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

m<-plot_muller(sim)
t<-plot_tissue(sim)
m+t
ggsave("tissue/muller.pdf", plot=m,height = 12, width = 10, dpi = 300, units = 'in')
n_w <- n_h <- 15
ncells <- 0.9 * n_w * n_h

# Sampling ncells with random box sampling of boxes of size n_w x n_h
bboxes <- sim$search_samples(c("Clone 4" = ncells), n_w, n_h,n_samples=3)
counter <- 1
for (bbox in bboxes) {
  sample_name <-paste0("SPN01_Sample_",counter)
  plot <- plot +
    ggplot2::geom_rect(xmin = bbox$lower_corner[1],
                       xmax = bbox$upper_corner[1],
                       ymin = bbox$lower_corner[2],
                       ymax = bbox$upper_corner[2],
                       fill = NA, color = "black")
  sim$sample_cells(sample_name, bbox$lower_corner, bbox$upper_corner)
  counter <- counter+1
}

forest <- sim$get_samples_forest()
plot_forest(forest) %>%
  annotate_forest(forest,samples = T,MRCAs = T)
#forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/races/Use_Cases/SPN01/samples_forest.sff")
forest$save("data/samples_forest.sff")
ggsave("tissue/forest.pdf", height = 12, width = 10, dpi = 300, units = 'in')

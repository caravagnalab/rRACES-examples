# runned in Orfeo, on Epyc node on 5th November 2024 
rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(patchwork)

seed <- 1679
set.seed(seed)

setwd('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN02/process')

sim = SpatialSimulation("SPN02", seed = seed, save_snapshots = TRUE)
sim$history_delta <- 1
sim$death_activation_level <- 50

# adding first clone
sim$add_mutant(name = "Clone 1", 
               growth_rates = 0.1,
               death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1000)

t1 = plot_tissue(sim)

# adding second clone 
sim$add_mutant(name = "Clone 2",
               growth_rates = 0.3,
               death_rates = 0.01)
sim$update_rates("Clone 1", rates = c(growth = .033))
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size("Clone 2", 1000)
sim$update_rates("Clone 1", rates = c(growth = 0, death = 0.025))
sim$run_up_to_size("Clone 2", 2000)
sim$update_rates("Clone 1", rates = c(growth = 0, death = 0.05))
sim$run_up_to_size("Clone 2", 4000)
sim$update_rates("Clone 1", rates = c(growth = 0, death = 0.1))
sim$run_up_to_size("Clone 2", 5000)

# most of cells from clone 1 are now dead
muller_clone2 = plot_muller(sim)
t2 = plot_tissue(sim)

# adding third clone
clone3_born = sim$get_clock()
saveRDS(object = clone3_born, file = "clone3_clock.rds")

sim$add_mutant(name = "Clone 3",
               growth_rates = 0.7,
               death_rates = 0.01)
sim$update_rates("Clone 2", rates = c(growth = .1))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 5000)
plot_muller(sim)

# now let's slowly kill clone 2
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.025))
sim$run_up_to_size("Clone 3", 6000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.05))
sim$run_up_to_size("Clone 3", 15000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.1))
sim$run_up_to_size("Clone 3", 25000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.15))
sim$run_up_to_size("Clone 3", 30000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.2))
sim$run_up_to_size("Clone 3", 40000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.3))
sim$run_up_to_size("Clone 3", 50000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.45))
sim$run_up_to_size("Clone 3", 70000)

# only clone 3 remains 
muller_clone3 = plot_muller(sim)
tissue_clone3 = plot_tissue(sim)

# SAMPLING
# two samples, at the same time and in two different places

bboxA_lower_corner <- c(460, 290)
bboxA_upper_corner <- c(500, 330)

bboxB_lower_corner <- c(390, 520)
bboxB_upper_corner <- c(430, 560)

# sample 1
sim$sample_cells("SPN02_1.1", bboxA_lower_corner, bboxA_upper_corner)

# sample 2
sim$sample_cells("SPN02_1.2", bboxB_lower_corner, bboxB_upper_corner)

plot_tissue(sim)
ggsave("tissue_sampled.png", bg = "white")

# Get forest ####
sampled_phylogeny <- sim$get_samples_forest()
sampled_phylogeny$save("samples_forest.sff")

plot_phylogeny <- plot_forest(sampled_phylogeny) %>%
  annotate_forest(sampled_phylogeny,
                  samples = TRUE,
                  MRCAs = TRUE,
                  exposures = FALSE,
                  facet_signatures = FALSE,
                  drivers = FALSE,
                  add_driver_label = FALSE)
ggsave('phylogeny.png', plot_phylogeny, bg = 'white')

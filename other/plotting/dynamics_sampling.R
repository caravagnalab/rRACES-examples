rm(list = ls())
library(ProCESS)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(ggpubr)
source('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/ProCESS-examples/plotting/utils_dynamics.R')

dir <- getwd()
set.seed(12345)
sim <- SpatialSimulation(name = 'SPN01', seed = 12345)
tissue_size <- sim$get_tissue_size()
sim$history_delta <- 1
sim$death_activation_level <- 50

tissue_plots <- list()
state_plots <- list()

clone_rates <- c(0.08,0.3,1,3)
names(clone_rates) <- c("Clone 1","Clone 2","Clone 3","Clone 4")

time_intervals <- 50


# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = 0.08, death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 500)


# Clone 2
sim$add_mutant(name = "Clone 2", growth_rates = 0.3, death_rates = 0.01)
sim$update_rates(species = "Clone 1", rates = c(growth = 0.06, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size(species = 'Clone 2', 700)

t1 <- get_timepoint_info(sim, 'T1')

sim$update_rates(species = "Clone 1", rates = c(growth = 0.02, death = 0.03))
sim$update_rates(species = "Clone 2", rates = c(growth = 0.4, death = 0.005))
sim$run_up_to_size(species = 'Clone 2', 1000)



# Grow 2
sim$update_rates(species = "Clone 1", rates = c(growth = 0.001, death = 0.3))
sim$run_up_to_size(species = 'Clone 2', 2000)
t2 <- get_timepoint_info(sim, 'T2')


# Clone 3
sim$add_mutant(name = "Clone 3", growth_rates = 1, death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.1, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 6000)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.05, death = 0.5))
sim$run_up_to_size("Clone 3", 8000)
t3 <- get_timepoint_info(sim, 'T3')

sim$update_rates(species = "Clone 2", rates = c(growth = 0.01, death = 0.8))
sim$run_up_to_size("Clone 3", 15000) ## orignal 15000
sim$update_rates(species = "Clone 2", rates = c(growth = 0.01, death = 1))
sim$run_up_to_size("Clone 3", 80000)

# Clone 4
sim$add_mutant(name = "Clone 4", growth_rates = 3 , death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.005, death = 1))
sim$update_rates(species = "Clone 3", rates = c(growth = 0.5, death = 0.02))
sim$mutate_progeny(sim$choose_cell_in("Clone 3"), "Clone 4")
sim$run_up_to_size("Clone 4", 10000) ## original 6000
sim$update_rates(species = "Clone 3", rates = c(growth = 0.02, death = 0.02))
sim$run_up_to_size("Clone 4", 60000) ## orignal 8000
t4 <- get_timepoint_info(sim, 'T4')

muller <- plot_muller(sim)
print("End simulation")

### Sampling###
bboxC <- sim$search_sample(c("Clone 4"= 1000),50, 50)
sim$sample_cells("SPN01_1.3", bboxC$lower_corner, bboxC$upper_corner)
s1.3 <- get_sampling_info(sim, "SPN01_1.3", bboxC)


bboxA <- sim$search_sample(c("Clone 4" = 500, 'Clone 3' = 500), 50,50)
sim$sample_cells("SPN01_1.1", bboxA$lower_corner, bboxA$upper_corner)
s1.1 <- get_sampling_info(sim, "SPN01_1.1", bboxA)


bboxB <- sim$search_sample(c("Clone 3" = 1000),50, 50)
sim$sample_cells("SPN01_1.2", bboxB$lower_corner, bboxB$upper_corner)
s1.2 <- get_sampling_info(sim, "SPN01_1.2", bboxB)


# Forest
forest <- sim$get_samples_forest()
forest$save("data/samples_forest.sff")

plt_forest <- plot_forest(forest) +
  theme(legend.position = "none")
piechart <- plot_state(sim) +
  theme(legend.position = "none")
timeseries <- plot_timeseries(sim) +
  theme(legend.position = "none")
muller <- plot_muller(sim) +
  theme(legend.position = "none")


# plot dynamics
sim_end <- sim$get_clock()

timepoints <- list(t1, t2, t3, t4)
times <- lapply(timepoints, FUN = function(t){
  t$time
}) %>% unlist()

labels <- lapply(timepoints, FUN = function(t){
  t$label
}) %>% unlist()

time.tb <- data.frame(Generation = times,
                      what =  labels,
                      event.type = "Timepoint")

time_plot <- ggplot(time.tb, aes(x = Generation, y = event.type, label = what)) +
  geom_line() +
  geom_point() +
  geom_text(hjust = 1.3, angle = 45) +
  xlim(0, sim_end) + 
  scale_y_discrete(name = "") +
  theme_minimal()

tissues <- lapply(timepoints, FUN = function(t){
  t$tissue
})

mullers <-  lapply(timepoints, FUN = function(t){
  t$muller
})

fires <- lapply(timepoints, FUN = function(t){
  t$fire
})

timeserie <-  lapply(timepoints, FUN = function(t){
  t$timeserie
})
states <- lapply(timepoints, FUN = function(t){
  t$state
})

legend_pos = 'None'
plt_tissues <-  patchwork::wrap_plots(tissues, nrow = 1, ncol = 4, guides = 'collect')  & theme(legend.position = legend_pos)
plt_states <-  patchwork::wrap_plots(states, nrow = 1, ncol = 4, guides = 'collect')  & theme(legend.position = legend_pos)
plt_mullers <-  patchwork::wrap_plots(tissues, nrow = 1, ncol = 4, guides = 'collect')  & theme(legend.position = legend_pos)
plt_fires <-   patchwork::wrap_plots(fires, nrow = 1, ncol = 4, guides = 'collect')  & theme(legend.position = legend_pos)
plt_timeseries <- patchwork::wrap_plots(timeserie, nrow = 1, ncol = 4, guides = 'collect')  & theme(legend.position = legend_pos)

v1 <- patchwork::wrap_plots(time_plot, plt_tissues, plt_states, plt_mullers, plt_timeseries, plt_fires, nrow = 6)
v2 <- patchwork::wrap_plots(time_plot, muller, timeseries + geom_vline(aes(xintercept = ifelse(.data$time %in% times, time, 0)), col = 'gray70'), plt_tissues, plt_states, nrow = 5)
v3 <- patchwork::wrap_plots(time_plot, plt_tissues, plt_states, plt_mullers, plt_timeseries, nrow = 5)

ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/dynamic_v1.pdf', v1, width = 210, height = 297, units = "mm")
ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/dynamic_v2.pdf', v2, width = 210, height = 297, units = "mm")
ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/dynamic_v3.pdf', v3, width = 210, height = 297, units = "mm")


# plot_sampling
sampling <- list(s1.1, s1.2, s1.3)
s_times <- lapply(sampling, FUN = function(t){
  t$time
}) %>% unlist()

s_labels <- lapply(sampling, FUN = function(t){
  t$label
}) %>% unlist()

sampling.tb <- data.frame(Generation = s_times,
                          what =  s_labels,
                          event.type = "Timepoint")

sample_time_plot <- ggplot(sampling.tb, aes(x = Generation, y = event.type, label = what)) +
  geom_line() +
  geom_point() +
  geom_text_repel(direction = "y",
                  point.padding = 0.5,
                  hjust = 0,
                  box.padding = 1,
                  seed = 123)  +
  scale_y_discrete(name = "") +
  theme_minimal() + 
  xlim(0, sim_end) 

samples_plot <- lapply(sampling, FUN = function(t){
  t$plot
})

legend_pos = 'None'
plt_samples <- patchwork::wrap_plots(samples_plot, nrow = 1, ncol = 4, guides = 'collect') & theme(legend.position = legend_pos)

clone_info <- get_clone_info(sim, forest)
plt_proportions <- patchwork::wrap_plots(plot_clonal_prop(clone_info, sim), nrow = 1, ncol = 4, guides = 'collect') & theme(legend.position = legend_pos)

final_sampling <- patchwork::wrap_plots(sample_time_plot, muller, plt_samples, plt_proportions, nrow = 4)
ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/sampling.pdf', final_sampling, width = 210, height = 297, units = "mm")

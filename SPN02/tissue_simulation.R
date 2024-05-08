library(rRACES)
library(tidyverse)
library(patchwork)

seed <- 1679
set.seed(seed)

sim = new(Simulation, 
          seed = seed,
          save_snapshot = F)
sim$update_tissue(1e3, 1e3)
sim$duplicate_internal_cells <- TRUE
sim$death_activation_level <- 50
dir.create("SPN02/tissue", recursive = T)
# adding first clone

sim$add_mutant(name = "Clone 1", 
               growth_rates = 0.5,
               death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1500)
clone2_born = sim$get_clock()

# adding second clone 

sim$add_mutant(name = "Clone 2",
               growth_rates = 1,
               death_rates = 0.01)

sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$update_rates("Clone 1", rates = c(growth = 0, death = 1))
sim$run_up_to_size("Clone 2", 3000)

plot_muller(sim)
ggsave("SPN02/tissue/clone2_muller.pdf")
plot_tissue(sim)
ggsave("SPN02/tissue/clone2_tissue.pdf")

# adding third clone

clone3_born = sim$get_clock()
sim$add_mutant(name = "Clone 3",
               growth_rates = 5,
               death_rates = 0.01)

sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
# sim$update_rates("Clone 1", rates = c(growth = 0, death = 3))
sim$update_rates("Clone 2", rates = c(growth = 0, death = 3.5))
sim$run_up_to_size("Clone 3", 15000)
plot_muller(sim)
ggsave(filename = "SPN02/tissue/clone3_muller.pdf")
plot_tissue(sim)
ggsave("SPN02/tissue/clone3_tissue.pdf")

# sampling twice from population of clone 3

box_w = 20
box1_p <- c(460, 420)
box1_q = box1_p + box_w

box2_p <- c(520, 430)
box2_q = box2_p + box_w

plot_tissue(sim) +
  geom_rect(xmin = box1_p[1], xmax = box1_q[1], 
            ymin = box1_p[2], ymax = box1_q[2], 
            fill = NA, color = "black") +
  geom_rect(xmin = box2_p[1], xmax = box2_q[1], 
            ymin = box2_p[2], ymax = box2_q[2], 
            fill = NA, color = "black") 

ggsave("SPN02/tissue/clone3_sampling.pdf")

sim$sample_cells("A", box1_p, box1_q)
sim$sample_cells("B", box2_p, box2_q)

# Get forest ####
sampled_phylogeny <- sim$get_samples_forest()
sampled_phylogeny$save("SPN02/tissue/forest.sff")

# Recap Plot ####
final_sampled_tissue <-plot_tissue(sim)
muller_plot <- plot_muller(sim)
piechart <- plot_state(sim)
timeseries_plot <- plot_timeseries(sim)
plot_phylogeny = plot_forest(sampled_phylogeny) %>%
  annotate_forest(sampled_phylogeny,
                  samples = TRUE,
                  MRCAs = TRUE,
                  exposures = FALSE,
                  facet_signatures = FALSE,
                  drivers = FALSE,
                  add_driver_label = FALSE) +
  ylim(0,10)

layout <- "
ACBD
#EE#
"
summary_patchwork = patchwork::wrap_plots(
  piechart, muller_plot, final_sampled_tissue, timeseries_plot,plot_phylogeny,
  guides = 'auto', design = layout
)

ggsave(filename = "SPN02/tissue/recap_plot.pdf", summary_patchwork, width = 15)


library(rRACES)
library(tidyverse)
library(patchwork)

seed <- 1679
set.seed(seed)
dt <- .1

sim = new(Simulation, 
          seed = seed,
          save_snapshot = T)
sim$update_tissue(1e3, 1e3)
sim$duplicate_internal_cells <- TRUE
sim$death_activation_level <- 50

# adding first clone

sim$add_mutant(name = "Clone 1", 
               growth_rates = 0.5,
               death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1500)
clone2_born = sim$get_clock()

# clone1_tissue = plot_tissue(sim)
# clone1_muller = plot_muller(sim)

# ggsave(clone1_tissue, file = "clone1_tissue.pdf")
# ggsave(clone1_muller, file = "clone1_muller_v2.pdf")

# update the growth rate of clone 1

# sim$update_rates("Clone 1", rates = c(growth = 0.7, death = 0.01))
# sim$run_up_to_size("Clone 1", 2000)

#sim$run_up_to_size("Clone 1", 1700)

# adding second clone 

sim$add_mutant(name = "Clone 2",
               growth_rates = 1,
               death_rates = 0.01)

sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$update_rates("Clone 1", rates = c(growth = 0, death = 0.2))
sim$run_up_to_size("Clone 2", 3000)

# sim$update_rates("Clone 1", rates = c(growth = 0, death = 0.02))
# sim$run_up_to_size("Clone 2", 3000)
plot_muller(sim)

# clone2_muller = plot_muller(sim)
# ggsave(clone2_muller, file = "clone2_muller_v1.pdf")
#sim$update_rates("Clone 2",rates = c(growth = 1.5, death=0.01))
#sim$update_rates("Clone 1",rates = c(growth = 0, death=0.5))
#sim$run_up_to_size("Clone 2",3000)

#sim$update_rates("Clone 1",rates = c(growth = 0, death=0.5))
#sim$run_up_to_size("Clone 2",5000)

#clone2_tissue <-plot_tissue(sim)
#clone2_muller <- plot_muller(sim)

#ggsave(clone2_tissue, file = "clone2_tissue.pdf")
#ggsave(clone2_muller, file = "clone2_muller_v2.pdf")

#sim$run_up_to_size("Clone 2",6000)

# adding third clone
clone3_born = sim$get_clock()
sim$add_mutant(name = "Clone 3",
               growth_rates = 5,
               death_rates = 0.01)

sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
#sim$run_up_to_size("Clone 3", 200)

#sim$update_rates("Clone 1", rates = c(growth = 0, death = 1))
# sim$update_rates("Clone 2", rates = c(growth = 0.7, death = 0.01))
# sim$run_up_to_size("Clone 2", 2000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.5))
#sim$update_rates("Clone 3", rates = c(growth = 5, death = 0.01))
sim$run_up_to_size("Clone 3", 8000)
# sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.2))

# sim$run_up_to_size("Clone 3", 15000)

plot_muller(sim)
# clone3_tissue <-plot_tissue(sim)
# ggsave(clone3_muller, file = "clone3_muller.pdf")
# ggsave(clone3_tissue, file = "clone3_tissue.pdf")

# killing the other clones

# sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.5))
# sim$update_rates("Clone 1", rates = c(growth = 0, death = 0.5))
# sim$run_up_to_size("Clone 3", 15000)

# sampling twice from population of clone 3

# box_w = box_h = 30
# ncells = 0.5*box_h*box_w
# bbox_1 <- sim$search_sample(c("Clone 3" = ncells), box_w, box_h)
# 
box_w = 20
# ncells = 0.8*box_h*box_w
# bbox_2 <- sim$search_sample(c("Clone 3" = ncells), 510, 635)
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

sim$sample_cells("A", box1_p, box1_q)
sim$sample_cells("B", box2_p, box2_q)

# Get forest ####
sampled_phylogeny <- sim$get_samples_forest()
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
  ylim(0,5)
layout <- "
ACBD
#EE#
"
summary_patchwork = patchwork::wrap_plots(
  piechart, muller_plot, final_sampled_tissue, timeseries_plot,plot_phylogeny,
  guides = 'auto', design = layout
)

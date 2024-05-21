library(rRACES)
# library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

seed <- 1679
set.seed(seed)

sim = new(Simulation, 
          seed = seed,
          save_snapshot = F)
sim$update_tissue(1e3, 1e3)
# sim$duplicate_internal_cells <- TRUE
sim$death_activation_level <- 50
dir.create("test/SPN02", recursive = T)
# adding first clone

sim$add_mutant(name = "Clone 1", 
               growth_rates = 0.5,
               death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1500)
# clone2_born = sim$get_clock()

# adding second clone 

sim$add_mutant(name = "Clone 2",
               growth_rates = 1,
               death_rates = 0.01)

sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$update_rates("Clone 1", rates = c(growth = 0, death = 1))
sim$run_up_to_size("Clone 2", 3000)

plot_muller(sim)
ggsave("test/SPN02/clone2_muller.pdf")
plot_tissue(sim)
ggsave("test/SPN02/clone2_tissue.pdf")

# adding third clone

clone3_born = sim$get_clock()
saveRDS(object = clone3_born, file = "test/SPN02/clone3_clock.rds")

sim$add_mutant(name = "Clone 3",
               growth_rates = 5,
               death_rates = 0.01)

sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
# sim$update_rates("Clone 1", rates = c(growth = 0, death = 3))
sim$update_rates("Clone 2", rates = c(growth = 0, death = 3.5))
sim$run_up_to_size("Clone 3", 15000)
plot_muller(sim)
ggsave(filename = "test/SPN02/clone3_muller.pdf")
plot_tissue(sim)
ggsave("test/SPN02/clone3_tissue.pdf")

# sampling twice from population of clone 3

# n_w <- n_h <- 30
# ncells <- 0.8 * n_w * n_h

# boxes <- sim$search_samples(c("Clone 3" = ncells), n_h, n_w, n_samples = 3)
# plot_s = plot_tissue(sim)
# for (bbox in boxes[c(1,3)]) {
#   plot_s <- plot_s +
#     ggplot2::geom_rect(xmin = bbox$lower_corner[1],
#                        xmax = bbox$upper_corner[1],
#                        ymin = bbox$lower_corner[2],
#                        ymax = bbox$upper_corner[2],
#                        fill = NA, color = "black")
# }

# plot_s
# ggsave(filename = "test/SPN02/sampling.pdf", plot = plot_s)
# bbox_1 <- sim$search_sample(c("Clone 3" = ncells), n_w, n_h)
# sim$sample_cells("S1", bbox_1$lower_corner, bbox_1$upper_corner)

# bbox_2 <- sim$search_sample(c("Clone 3" = ncells), n_w, n_h)
# sim$sample_cells("S2", bbox_2$lower_corner, bbox_2$upper_corner)

box_w = 45
box1_p <- c(430, 430)
box1_q = box1_p + box_w

box2_p <- c(500, 420)
box2_q = box2_p + box_w

# box_1 = bbox_sampler(sim, which="Clone 3", n = 2000, n_w = 50, n_h = 50)
# box_2 = bbox_sampler(sim, which="Clone 3", n = 2000, n_w = 50, n_h = 50)

plot_tissue(sim) +
  geom_rect(xmin = box1_p[1], xmax = box1_q[1], 
            ymin = box1_p[2], ymax = box1_q[2], 
            fill = NA, color = "#C40C0C") +
  geom_rect(xmin = box2_p[1], xmax = box2_q[1], 
            ymin = box2_p[2], ymax = box2_q[2], 
            fill = NA, color = "#333A73")

ggsave("test/SPN02/sampling.pdf")

sim$sample_cells("A", box1_p, box1_q)
sim$sample_cells("B", box2_p, box2_q)
plot_tissue(sim)
ggsave("test/SPN02/tissue_sampled.pdf")
# Get forest ####
sampled_phylogeny <- sim$get_samples_forest()
sampled_phylogeny$save("test/SPN02/samples_forest.sff")

# Recap Plot ####
final_sampled_tissue <-plot_tissue(sim)
muller_plot <- plot_muller(sim)
piechart <- plot_state(sim)
timeseries_plot <- plot_timeseries(sim)
plot_phylogeny <- plot_forest(sampled_phylogeny) %>%
  annotate_forest(sampled_phylogeny,
                  samples = TRUE,
                  MRCAs = TRUE,
                  exposures = FALSE,
                  facet_signatures = FALSE,
                  drivers = FALSE,
                  add_driver_label = FALSE) +
  ylim(0,7)
print("created phylogeny")

plot_phylogeny
# ggsave(object = plot_phylogeny, filename = "test/SPN02/phylogeny.pdf", width = 10)

print("saved plot phylo")

layout <- "
#AA#BB#
#CC#DD#
#EEEEE#
"
summary_patchwork = patchwork::wrap_plots(
  piechart, muller_plot, final_sampled_tissue, timeseries_plot,plot_phylogeny,
  guides = 'auto', design = layout
)

ggsave(filename = "test/SPN02/recap_plot2.pdf", summary_patchwork, width = 15, height = 15)
# runned in Orfeo, on Epyc node on 5th November 2024 
rm(list=ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
library(english)
source("rRACES-examples/SPN02/plotting_sampling.R")

setwd("SPN02/results/")

seed <- 1679
set.seed(seed)

sim = SpatialSimulation("SPN02")
sim$history_delta <- 1
# sim$border_growth_model = TRUE
# sim$duplicate_internal_cells <- TRUE
sim$death_activation_level <- 50

# adding first clone
sim$add_mutant(name = "Clone 1", 
               growth_rates = 0.1,
               death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1000)

t1 = plot_tissue(sim)
ggsave("tissue_clone1.png", plot = t1)

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
sim$update_rates("Clone 1", rates = c(growth = 0, death = 0.08))
sim$run_up_to_size("Clone 2", 5000)

# most of cells from clone 1 are now dead
muller_clone2 = plot_muller(sim)
ggsave("muller_clone2.png", plot = muller_clone2)

t2 = plot_tissue(sim)
ggsave("tissue_clone2.png", plot = t2)

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
ggsave(filename = "clone3_muller.png")

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
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.25))
sim$run_up_to_size("Clone 3", 50000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = 0.35))
sim$run_up_to_size("Clone 3", 65000)

# only clone 3 remains 
muller_clone3 = plot_muller(sim)
ggsave(filename = "clone3_muller.png", plot = muller_clone3)
tissue_clone3 = plot_tissue(sim)
ggsave("clone3_tissue.png", plot = tissue_clone3)

# sample the tissue
# 1. setting the boxes --> doing this is equivalent of doing the stuff below, but comparable with the automatic search boxes 
bbox = new(TissueRectangle, c(400, 350), 35, 35)
bbox2 = new(TissueRectangle, c(350, 530), 35, 35)

# plotting_sample = function(sim, samples_timing, boxes) {
#   
#   color_type = c("Polyclonal" = "#C40C0C", "Clonal" = "#333A73")
#   
#   layout = "
#     AABBCCDD
#     "
#   
#   p_list <- lapply(english::ordinal(seq(1:4)), function(x) {ggplot() + ggtitle(paste(x, "sampling"))})
#   
#   for (t in seq_along(samples_timing)) {
#     samples = samples_timing[[t]]
#     print(samples)
#     p_list[[t]] <- rRACES::plot_tissue(sim)
#     for (i in seq_along(samples)) {
#       
#       sample = samples[i]
#       bbox = boxes[[sample]]
#       
#       # print(sample)
#       sample_data = tibble(type = names(sample), 
#                            sample = unname(sample), 
#                            x_min = bbox$lower_corner[1],
#                            x_max = bbox$upper_corner[1],
#                            y_min = bbox$lower_corner[2],
#                            y_max = bbox$upper_corner[2])
#       
#       p_list[[t]] =  p_list[[t]] +
#         geom_rect(data = sample_data, aes(color = type, 
#                                           xmin = x_min,
#                                           xmax = x_max,
#                                           ymin = y_min,
#                                           ymax = y_max),
#                   fill = NA, inherit.aes = FALSE) +
#         scale_color_manual(values = color_type) +
#         ggtitle(paste(english::ordinal(t), "sampling")) +
#         guides(color = guide_legend(nrow = 2, title = "Type of sample", title.position = "top"))
#     }
#   }
#   
#   p = wrap_plots(p_list, design = layout) + plot_layout(guides = 'collect')
#   return(p)
# }

s = c("Clonal" = "SPN02_1.1", "Clonal" = "SPN02_1.2")
timing = list("T1" = s)
box = list("SPN02_1.1" = bbox, "SPN02_1.2" = bbox2)
sampling_plot = plotting_sample(sim, samples_timing = timing, boxes = box)
ggsave("sampling.png", bg = "white", plot = sampling_plot, width = 12, height = 5)

lapply(box %>% names, function(x) {
  sim$sample_cells(x, box[[x]]$lower_corner, box[[x]]$upper_corner)
})
##############################

# box_w = 20
# box1_p <- c(430, 350)
# box1_q = box1_p + box_w
# 
# box2_p <- c(300, 510)
# box2_q = box2_p + box_w

# box_1 = bbox_sampler(sim, which="Clone 3", n = 2000, n_w = 50, n_h = 50)
# box_2 = bbox_sampler(sim, which="Clone 3", n = 2000, n_w = 50, n_h = 50)



# sampling_plot = plot_tissue(sim) +
#   geom_rect(xmin = box1_p[1], xmax = box1_q[1], 
#             ymin = box1_p[2], ymax = box1_q[2], 
#             fill = NA, color = "#C40C0C") +
#   geom_rect(xmin = box2_p[1], xmax = box2_q[1], 
#             ymin = box2_p[2], ymax = box2_q[2], 
#             fill = NA, color = "#333A73")

# sim$sample_cells("A", box1_p, box1_q)
# sim$sample_cells("B", box2_p, box2_q)
plot_tissue(sim)
ggsave("tissue_sampled.png", bg = "white")
# Get forest ####
sampled_phylogeny <- sim$get_samples_forest()
sampled_phylogeny$save("samples_forest.sff")

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
                  add_driver_label = FALSE) #+
  # ylim(0,7)
print("created phylogeny")

# plot_phylogeny
# ggsave(object = plot_phylogeny, filename = "test/SPN02/phylogeny.pdf", width = 10)

print("saved plot phylo")

layout <- "
AA#BB#CC
DDDDEEEE
FF#GG#HH
IIIIIIII
IIIIIIII
"
summary_patchwork = patchwork::wrap_plots(
  t1, t2, tissue_clone3,
  muller_clone2, muller_plot, 
  piechart,timeseries_plot, sampling_plot,
  plot_phylogeny,
  guides = 'auto', design = layout
)

ggsave(filename = "summary.png", summary_patchwork, width = 15, height = 20)

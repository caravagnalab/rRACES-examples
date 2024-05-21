
# Recap Plot ####
muller_plot <- plot_muller(sim)
piechart <- plot_state(sim)
timeseries_plot <- plot_timeseries(sim)
final_sampled_tissue <- plot_tissue(sim)

#plot_phylogeny = plot_forest(sampled_phylogeny) %>% 
#  annotate_forest(
#    sampled_phylogeny, 
#    samples = TRUE, 
#    MRCAs = TRUE, 
#    exposures = FALSE, 
#    facet_signatures = FALSE, 
#    drivers = FALSE, 
#    add_driver_label = FALSE
#    )

layout <- "
ACBD
#EE#
"

summary_patchwork = patchwork::wrap_plots(
  piechart, 
  muller_plot, 
  final_sampled_tissue, 
  timeseries_plot, 
  #plot_phylogeny, 
  guides = 'auto', 
  design = layout
)

summary_patchwork



#===============================================================================
# STORAGE
#===============================================================================

# sample A
bbox = tibble(lower_corner = c(1000,1000), upper_corner = c(1040,1080))

plot_tissue(sim, num_of_bins  = 50) + 
  geom_rect(
    xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
    ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
    fill = NA, color = "black"
  )

sim$sample_cells("AA", bbox$lower_corner, bbox$upper_corner)
forest <- sim$get_samples_forest()
plot_forest(forest)
forest$get_nodes() %>% filter(!is.na(sample)) %>% group_by(mutant) %>% summarize(n = length(cell_id))


# sample B
bbox = tibble(lower_corner = c(1300,1060), upper_corner = c(1400,1080))

plot_tissue(sim, num_of_bins  = 50) + 
  geom_rect(
    xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
    ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
    fill = NA, color = "black"
  )

sim$sample_cells("B", bbox$lower_corner, bbox$upper_corner)
forest <- sim$get_samples_forest()
plot_forest(forest)
forest$get_nodes() %>% filter(!is.na(sample), sample == "B") %>% group_by(mutant) %>% summarize(n = length(cell_id))


bbox = tibble(lower_corner = c(1000,1000), upper_corner = c(1040,1080))

plot_tissue(sim, num_of_bins  = 50) + 
  geom_rect(
    xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
    ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
    fill = NA, color = "black"
  )

sim$sample_cells("C", bbox$lower_corner, bbox$upper_corner)
forest <- sim$get_samples_forest()
plot_forest(forest)
forest$get_nodes() %>% filter(!is.na(sample)) %>% group_by(mutant) %>% summarize(n = length(cell_id))

# sample D
bbox = tibble(lower_corner = c(1000,1000), upper_corner = c(1040,1080))

plot_tissue(sim, num_of_bins  = 50) + 
  geom_rect(
    xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
    ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
    fill = NA, color = "black"
  )

sim$sample_cells("D", bbox$lower_corner, bbox$upper_corner)
forest <- sim$get_samples_forest()
plot_forest(forest)
forest$get_nodes() %>% filter(!is.na(sample)) %>% group_by(mutant) %>% summarize(n = length(cell_id))


# sample E
bbox = tibble(lower_corner = c(1000,1000), upper_corner = c(1040,1080))

plot_tissue(sim, num_of_bins  = 50) + 
  geom_rect(
    xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
    ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
    fill = NA, color = "black"
  )

sim$sample_cells("E", bbox$lower_corner, bbox$upper_corner)
forest <- sim$get_samples_forest()
plot_forest(forest)
forest$get_nodes() %>% filter(!is.na(sample)) %>% group_by(mutant) %>% summarize(n = length(cell_id))



sim$run_up_to_time(sim$get_clock() + 2)

sim$run_up_to_size("C2", 10)

sim$run_up_to_event("death", "C2", 100)

sim$run_up_to_time(sim$get_clock() + 1)




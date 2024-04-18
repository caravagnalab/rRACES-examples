
# ==============================================================================
# =================================== SPN06 ====================================
# ==============================================================================

# To-Do
# - check the coordinates

#devtools::install_github("caravagnalab/rRACES")
#devtools::load_all("/Users/azadsadr/Documents/packages/rRACES")

# loading the library
library(dplyr)
library(ggplot2)
library(rRACES)


#-------------------------------------------------------------------------------
#---------------------------- create the tissue --------------------------------
#-------------------------------------------------------------------------------
sim <- new(Simulation, "SPN06", seed = 3, save_snapshot = F)

sim$duplicate_internal_cells <- T

sim$update_tissue("Lung", 2e3, 2e3)


#-------------------------------------------------------------------------------
#------------------------- Rising of Clones: [C1, C2, C3, C4] ------------------
#-------------------------------------------------------------------------------

# clone 01
sim$add_mutant(name = "C1", growth_rates = 2, death_rates = 0)
sim$place_cell("C1", 1000, 1000)
sim$run_up_to_size("C1", 1e2)

# clone 02
sim$add_mutant(name = "C2", growth_rates = 4, death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("C1"), "C2")
sim$run_up_to_size("C2", 1e2)


# clone 03
sim$add_mutant(name = "C3", growth_rates = 6, death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C3")
sim$run_up_to_size("C3",1e2)

# clone 04
sim$add_mutant(name = "C4", growth_rates = 8, death_rates = 0)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C4")
sim$run_up_to_size("C4",1e5)

print(sim)
plot_state(sim)
plot_tissue(sim)
plot_tissue(sim,num_of_bins  = 50)
plot_tissue(sim, num_of_bins = 250) + facet_wrap(~species)
plot_timeseries(sim)
plot_muller(sim)

#-------------------------------------------------------------------------------
#------------------------------ Sampling [A, B] --------------------------------
#-------------------------------------------------------------------------------

# ======= START TEST

bbox = sim$search_sample(c("C3" = 100), nw = 50, nh = 50)
sim$sample_cells("A", bbox$lower_corner, bbox$upper_corner)
forest <- sim$get_samples_forest()

bbox = sim$search_sample(c("C2" = 100,"C3" = 100), nw = 50, nh = 50)
sim$sample_cells("B", bbox$lower_corner, bbox$upper_corner)
forest <- sim$get_samples_forest()

# ======= END TEST

# sample A
bbox = tibble(lower_corner = c(1000,1000), upper_corner = c(1040,1080))

plot_tissue(sim, num_of_bins  = 50) + 
  geom_rect(
    xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
    ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
    fill = NA, color = "black"
    )

sim$sample_cells("A", bbox$lower_corner, bbox$upper_corner)
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



#-------------------------------------------------------------------------------
#------------------------------ Treatment ENTERS I -----------------------------
#-------------------------------------------------------------------------------

sim$update_rates("C1", c(death = 125))
sim$update_rates("C2", c(death = 125))
sim$update_rates("C3", c(death = 125))
sim$update_rates("C4", c(death = 125))
sim$run_up_to_time(sim$get_clock() + 1)

#-------------------------------------------------------------------------------
#----------------------------- Rising of Clones: [C5] ----------------------------
#-------------------------------------------------------------------------------

sim$add_mutant(name = "C5", growth_rates = 120, death_rates = 0)
sim$mutate_progeny(
  as.data.frame((sim$get_cells() %>% filter(mutant == 'C4') %>% arrange(position_x))[1,]), 
  "C5"
  )
sim$run_up_to_size("C5", 1e4)
plot_tissue(sim)

#-------------------------------------------------------------------------------
#-------------------------------- Sampling [C] ---------------------------------
#-------------------------------------------------------------------------------

# START TEST
#n_w <- n_h <- 15
#ncells <- 1*n_w*n_h
#bbox <- sim$search_sample(c("C5" = ncells), n_w, n_h)
#sim$sample_cells("C", bbox$lower_corner, bbox$upper_corner)
# END TEST


# sample C
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


#-------------------------------------------------------------------------------
#----------------------------- Treatment ENTERS II -----------------------------
#-------------------------------------------------------------------------------

sim$update_rates("C5", c(death = 125))
sim$run_up_to_time(sim$get_clock() + 1)


#-------------------------------------------------------------------------------
#----------------------------- Rising of Clones [C6]  ---------------------------
#-------------------------------------------------------------------------------

sim$add_mutant(name = "C6", growth_rates = 120, death_rates = 0)
sim$mutate_progeny(
  as.data.frame((sim$get_cells() %>% filter(mutant == 'C5') %>% arrange(position_x))[1,]), 
  "C6"
)
sim$run_up_to_size("C6", 1e4)
plot_tissue(sim)


#-------------------------------------------------------------------------------
#----------------------------- Sampling [D, E] ---------------------------------
#-------------------------------------------------------------------------------

# START TEST
#n_w <- n_h <- 15
#ncells <- 1 * n_w * n_h
#bbox <- sim$search_sample(c("C6" = ncells), n_w, n_h)
#sim$sample_cells("D", bbox$lower_corner, bbox$upper_corner)

#n_w <- n_h <- 15
#ncells <- 1*n_w*n_h
#bbox <- sim$search_sample(c("C6" = ncells), n_w, n_h)
#sim$sample_cells("E", bbox$lower_corner, bbox$upper_corner)
# END TEST


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




#-------------------------------------------------------------------------------


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









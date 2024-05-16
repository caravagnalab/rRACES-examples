
# ==============================================================================
# =================================== SPN06 ====================================
# ==============================================================================

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
sim$update_tissue("Lung", 1000, 1000)

#-------------------------------------------------------------------------------
#---------------------- Rising of Clones I : [C1, C2, C3, C4] ------------------
#-------------------------------------------------------------------------------

# clone 01
sim$add_mutant(name = "C1", growth_rates = 0.1, death_rates = 0.01)
sim$place_cell("C1", 500, 500)
sim$run_up_to_size("C1", 1000)

#plot_tissue(sim)
#sim$get_counts()

# clone 02
sim$add_mutant(name = "C2", growth_rates = 0.5, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C1"), "C2")
sim$run_up_to_size("C2", 1000)

#plot_tissue(sim)
#sim$get_counts()


# clone 03
sim$add_mutant(name = "C3", growth_rates = 1, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C3")
sim$run_up_to_size("C3",1000)

#plot_tissue(sim)
#sim$get_counts()


# clone 04
sim$add_mutant(name = "C4", growth_rates = 3, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C4")
sim$run_up_to_size("C4",8000)

#plot_tissue(sim)
#sim$get_counts()

#-------------------------------------------------------------------------------
#------------------------------ Sampling [A, B] --------------------------------
#-------------------------------------------------------------------------------

# sample A
bbox = sim$search_sample(c("C2" = 20, "C3" = 180), nw = 20, nh = 20)
sim$sample_cells("A", bbox$lower_corner, bbox$upper_corner)

# sample B
bbox = sim$search_sample(c("C2" = 120, "C3" = 80), nw = 20, nh = 20)
sim$sample_cells("B", bbox$lower_corner, bbox$upper_corner)


#-------------------------------------------------------------------------------
#------------------------------ Treatment ENTERS I -----------------------------
#-------------------------------------------------------------------------------

chemo1_start <- sim$get_clock()

sim$update_rates("C1", rates = c(growth = 0, death = 3))
sim$update_rates("C2", rates = c(growth = 0, death = 3))
sim$update_rates("C3", rates = c(growth = 0, death = 3))
sim$update_rates("C4", rates = c(growth = 0, death = 3))

while ((sim$get_cells() %>% dplyr::filter(mutant == "C2") %>% nrow()) > 1000) {
  sim$run_up_to_time(sim$get_clock() + 0.1)
}

chemo1_end <- sim$get_clock()

#-------------------------------------------------------------------------------
#--------------------------- Rising of Clones II : [C5] ------------------------
#-------------------------------------------------------------------------------

sim$add_mutant(name = "C5", growth_rates = 4, death_rates = 0.01)

sim$mutate_progeny(
  as.data.frame((sim$get_cells() %>% filter(mutant == 'C2') %>% arrange(position_x))[1,]), 
  "C5"
  )

sim$run_up_to_size("C5", 1e5)

#-------------------------------------------------------------------------------
#-------------------------------- Sampling [C] ---------------------------------
#-------------------------------------------------------------------------------

# sample C
bbox = sim$search_sample(c("C5" = 100), nw = 20, nh = 20)
sim$sample_cells("C", bbox$lower_corner, bbox$upper_corner)


#-------------------------------------------------------------------------------
#----------------------------- Treatment ENTERS II -----------------------------
#-------------------------------------------------------------------------------

chemo2_start <- sim$get_clock()

sim$update_rates("C5", rates = c(growth = 0, death = 3))

while ((sim$get_cells() %>% dplyr::filter(mutant == "C5") %>% nrow()) > 1000) {
  sim$run_up_to_time(sim$get_clock() + 0.1)
}

chemo2_end <- sim$get_clock()

#-------------------------------------------------------------------------------
#------------------------- Rising of Clones III : [C6]  -----------------------s
#-------------------------------------------------------------------------------

sim$add_mutant(name = "C6", growth_rates = 4, death_rates = 0)

sim$mutate_progeny(
  as.data.frame((sim$get_cells() %>% filter(mutant == 'C5') %>% arrange(position_x))[1,]), 
  "C6"
)
sim$run_up_to_size("C6", 1e5)


#-------------------------------------------------------------------------------
#----------------------------- Sampling [D, E] ---------------------------------
#-------------------------------------------------------------------------------

# sample C
bbox = sim$search_sample(c("C6" = 200), nw = 20, nh = 20)
sim$sample_cells("D", bbox$lower_corner, bbox$upper_corner)

# sample C
bbox = sim$search_sample(c("C6" = 200), nw = 20, nh = 20)
sim$sample_cells("E", bbox$lower_corner, bbox$upper_corner)


#-------------------------------------------------------------------------------
#----------------------------------- Save --------------------------------------
#-------------------------------------------------------------------------------

forest <- sim$get_samples_forest()
forest$save("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/forest.sff")

chemo_timing <- list(
  chemo1_start = round(chemo1_start, 0), 
  chemo1_end = round(chemo1_end, 0), 
  chemo2_start = round(chemo2_start, 0), 
  chemo2_end = round(chemo2_end, 0)
)
saveRDS(chemo_timing, "/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/chemo_timing.rds")


#-------------------------------------------------------------------------------
#----------------------------------- Plots -------------------------------------
#-------------------------------------------------------------------------------

print(sim)
sim$get_cells()
sim$get_counts()
forest$get_samples_info()
forest$get_nodes() %>% filter(!is.na(sample)) %>% group_by(mutant) %>% summarize(n = length(cell_id))

plot_state(sim)
plot_timeseries(sim)


#--- [MULLER PLOT] -------------------------------------------------------------
muller <- plot_muller(sim)
ggsave("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/muller_04.pdf", muller, dpi=300, width = 16, height = 8)

#--- [TISSUE PLOT] -------------------------------------------------------------
tissue <- plot_tissue(sim)
#plot_tissue(sim, num_of_bins  = 50)
#plot_tissue(sim, num_of_bins = 250) + facet_wrap(~species)
ggsave("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/tissue_04.pdf", tissue, dpi=300, width = 16, height = 8)

#--- [FOREST PLOT] -------------------------------------------------------------
plot_forest(forest) %>%
  annotate_forest(forest, 
                  samples = TRUE, 
                  MRCAs = TRUE, 
                  exposures = FALSE, 
                  facet_signatures = FALSE, 
                  drivers = FALSE, 
                  add_driver_label = FALSE)
ggsave("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/forest.pdf", dpi=300, width = 12, height = 8)












# ==============================================================================
# =================================== SPN06 ====================================
# ==============================================================================


#devtools::install_github("caravagnalab/rRACES")
#devtools::load_all("/u/cdslab/ahaghighi/scratch/packages/rRACES")
#devtools::load_all("~/Documents/packages/rRACES")

.libPaths("/orfeo/cephfs/scratch/cdslab/ahaghighi/Rlibs")

library(dplyr)
library(ggplot2)
library(rRACES)
library(patchwork)

seed <- 2024
set.seed(seed)


#-------------------------------------------------------------------------------
#---------------------------- create the tissue --------------------------------
#-------------------------------------------------------------------------------
sim <- SpatialSimulation(name = "SPN06", width = 2000, height = 2000, save_snapshots = FALSE, seed = seed)
sim$history_delta <- 1
sim$death_activation_level <- 50


#-------------------------------------------------------------------------------
#---------------------- Rising of Clones I : [C1, C2, C3, C4] ------------------
#-------------------------------------------------------------------------------

#--------------------------------- clone 01 ------------------------------------
sim$add_mutant(name = "C1", growth_rates = 0.1, death_rates = 0.01) # add species
sim$place_cell("C1", 1000, 1000) # displace the initial cell in the tissue
sim$run_up_to_size("C1", 1000)

#--------------------------------- clone 02 ------------------------------------
sim$add_mutant(name = "C2", growth_rates = 0.3, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C1"), "C2")

sim$update_rates("C1", rates = c(growth = 0.05, death = 0.01))
sim$run_up_to_size("C2", 1000)

sim$update_rates("C1", rates = c(growth = 0, death = 0.025))
sim$run_up_to_size("C2", 2000)

sim$update_rates("C1", rates = c(growth = 0, death = 0.05))
sim$update_rates("C2", rates = c(growth = 0.2, death = 0.02))
sim$run_up_to_size("C2", 8000)

#--------------------------------- clone 03 ------------------------------------
sim$add_mutant("C3",growth_rates = 0.5, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C3")
sim$update_rates("C2", rates = c(growth = 0.2, death = 0.06))
sim$run_up_to_size("C3", 1000)

#--------------------------------- clone 04 ------------------------------------
sim$add_mutant(name = "C4", growth_rates = 0.5, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C4")
sim$update_rates("C3", rates = c(growth = 0.2, death = 0.01))
sim$run_up_to_size("C4", 1000)


#-------------------------------------------------------------------------------
#------------------------------ Sampling [A, B] --------------------------------
#-------------------------------------------------------------------------------
bbox_width <- 20

# Box A
bbox1_p <- c(970, 895) # OK
bbox1_q <- bbox1_p + bbox_width

# Box B
bbox2_p <- c(900, 925) # 80% 0% 20%
bbox2_q <- bbox2_p + bbox_width

sim$sample_cells("A", bottom_left = bbox1_p, top_right = bbox1_q)
sim$sample_cells("B", bottom_left = bbox2_p, top_right = bbox2_q)
message("Samples A and B Extracted!")
sample_AB <- plot_tissue(sim)


#-------------------------------------------------------------------------------
#---------------------------- CHEMOTHERAPY ENTERS I ----------------------------
#-------------------------------------------------------------------------------
chemo1_start <- sim$get_clock()

sim$update_rates("C1", rates = c(growth = 0, death = 1))
sim$update_rates("C2", rates = c(growth = 0, death = 1))
sim$update_rates("C3", rates = c(growth = 0, death = 1))
sim$update_rates("C4", rates = c(growth = 0, death = 1))

#while ((sim$get_cells() %>% filter(mutant == "C2") %>% nrow()) > 100) {
#  sim$run_up_to_time(sim$get_clock() + 0.1)
#}
sim$run_until(sim$var("C2")<100)


chemo1_end <- sim$get_clock()


#-------------------------------------------------------------------------------
#--------------------------- Rising of Clones II : [C5] ------------------------
#-------------------------------------------------------------------------------
sim$add_mutant(name = "C5", growth_rates = 0.2, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C5")
sim$run_up_to_size("C5", 16000)


#-------------------------------------------------------------------------------
#-------------------------------- Sampling [C] ---------------------------------
#-------------------------------------------------------------------------------
bbox_c = sim$search_sample(c("C5" = 360), nw = 20, nh = 20)
sim$sample_cells("C", bbox_c$lower_corner, bbox_c$upper_corner)
message("Sample C Extracted!")
sample_C <- plot_tissue(sim)


#-------------------------------------------------------------------------------
#----------------------------- CHEMOTHERAPY ENTERS II --------------------------
#-------------------------------------------------------------------------------
chemo2_start <- sim$get_clock()

sim$update_rates("C5", rates = c(growth = 0, death = 3))

while ((sim$get_cells() %>% dplyr::filter(mutant == "C5") %>% nrow()) > 100) {
  sim$run_up_to_time(sim$get_clock() + 0.1)
}

chemo2_end <- sim$get_clock()


#-------------------------------------------------------------------------------
#------------------------ Rising of Clones III : [C6, C7] ----------------------
#-------------------------------------------------------------------------------

#--------------------------------- clone 06 ------------------------------------
sim$add_mutant(name = "C6", growth_rates = 0.2, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C5"), "C6")
sim$run_up_to_size("C6", 3000)
message("clone 06 created!")

#--------------------------------- clone 07 ------------------------------------
sim$add_mutant("C7", growth_rates = 0.2, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C6"), "C7")
sim$update_rates("C6", rates = c(growth = 0.1, death = 0.03))
sim$run_up_to_size("C7", 2000)


#-------------------------------------------------------------------------------
#----------------------------- Sampling [D, E] ---------------------------------
#-------------------------------------------------------------------------------
bbox_d = sim$search_sample(c("C6" = 360), nw = 20, nh = 20)
sim$sample_cells("D", bbox_d$lower_corner, bbox_d$upper_corner)

bbox_e = sim$search_sample(c("C7" = 360), nw = 20, nh = 20)
sim$sample_cells("E", bbox_e$lower_corner, bbox_e$upper_corner)

message("Samples D and E Extracted!")
sample_DE <- plot_tissue(sim)


#-------------------------------------------------------------------------------
#----------------------------------- PLOTS -------------------------------------
#-------------------------------------------------------------------------------
state <- plot_state(sim)
muller <- plot_muller(sim)
ts <- plot_timeseries(sim)
forest <- sim$get_samples_forest()
forest_plot <- plot_forest(forest) %>% annotate_forest(forest)

design <- "abcde
           fff##
           ggg##
           ggg##"

p <- patchwork::wrap_plots(
  a = sample_AB, 
  b = sample_C, 
  c = sample_DE, 
  d = state, 
  e = ts, 
  f = muller, 
  g = forest_plot,
  guides = "auto", 
  design = design
)

message("final plot created!")


#-------------------------------------------------------------------------------
#----------------------------------- Save --------------------------------------
#-------------------------------------------------------------------------------
setwd("/u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06")

curr_dir <- getwd()

print(paste0("curr_dir: ", curr_dir))

ggsave(filename = paste0(curr_dir, "/plots/final.png"), plot = p, width = 16, height = 10, dpi = 300)
#ggsave(filename = "~/scratch/packages/rRACES-examples/SPN06/plots/final2.png", plot = p, width = 16, height = 10, dpi = 300)
#ggsave(filename = paste(curr_dir, "/plots/final_test2.png", sep = ""), plot = p, width = 16, height = 10, dpi = 300)

forest$save( paste(curr_dir, "/data/forest.sff", sep = "") )

chemo_timing <- list(
  chemo1_start = chemo1_start, # round(chemo1_start, 0)
  chemo1_end = chemo1_end,     # round(chemo1_end, 0)
  chemo2_start = chemo2_start, # round(chemo2_start, 0)
  chemo2_end = chemo2_end      # round(chemo2_end, 0)
)

saveRDS( chemo_timing, paste(curr_dir, "/data/chemo_timing.rds", sep = "") )

message("DONE!, everything safe and sound")



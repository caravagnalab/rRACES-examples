
#-------------------------------------------------------------------------------
#------------------------------ Setting up environment -------------------------
#-------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) # clears all objects including hidden objects from the workspace
gc() # free-up memory (garbage collector)

#devtools::install_github("caravagnalab/ProCESS")
#path <- "/orfeo/cephfs/home/cdslab/ahaghighi/R/x86_64-pc-linux-gnu-library/4.5" # new library path

library(ProCESS)
library(dplyr)
library(ggplot2)
library(patchwork)


dir <- getwd()
outdir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN06/process/"
setwd(outdir)
seed = 19999
set.seed(seed)


# ==============================================================================
# =================================== SPN06 ====================================
# ==============================================================================


#-------------------------------------------------------------------------------
#---------------------------- create the tissue --------------------------------
#-------------------------------------------------------------------------------
sim <- SpatialSimulation(name = "SPN06", width = 1000, height = 1000, save_snapshots = FALSE, seed = seed)
sim$history_delta <- 1 # 1 or 0.1
sim$death_activation_level <- 50
#sim$border_growth_model <- FALSE ?????


#-------------------------------------------------------------------------------
#---------------------- Rising of Clones I : [C1, C2, C3, C4] ------------------
#-------------------------------------------------------------------------------

#--------------------------------- clone 01 ------------------------------------
sim$add_mutant(name = "C1", growth_rates = 0.1, death_rates = 0.01) # add species
sim$place_cell("C1", 500, 500) # displace the initial cell in the tissue
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
sim$update_rates("C2", rates = c(growth = 0.2, death = 0.08))
sim$run_up_to_size("C3", 1000)

#--------------------------------- clone 04 ------------------------------------
sim$add_mutant(name = "C4", growth_rates = 0.6, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C4")
sim$update_rates("C2", rates = c(growth = 0.2, death = 0.1))
sim$update_rates("C3", rates = c(growth = 0.2, death = 0.01))
sim$run_up_to_size("C4", 4000)


#-------------------------------------------------------------------------------
#------------------------------ Sampling [A, B] --------------------------------
#-------------------------------------------------------------------------------

#--------------------------------- sample 1.1 ----------------------------------
a_lowercorbner <- as.integer(c(405, 420))
a_uppercorner <- as.integer(c(454, 469))
sim$get_cells(a_lowercorbner, a_uppercorner) %>% group_by(mutant) %>% summarise(count = n()) 
# C2(156) 6%, C3(2318) 94%, total(2474)

SPN_1.1_plot <- plot_tissue(sim, num_of_bins=100) +
  ggplot2::geom_rect(ggplot2::aes(
    xmin=a_lowercorbner[1], 
    xmax=a_uppercorner[1], 
    ymin=a_lowercorbner[2], 
    ymax=a_uppercorner[2]), 
    color= 'red', fill='white')

sim$sample_cells("SPN06_1.1", a_lowercorbner, a_uppercorner)


#--------------------------------- sample 1.2 ----------------------------------
b_lowercorbner <- as.integer(c(440, 380))
b_uppercorner <- as.integer(c(489, 429))
sim$get_cells(b_lowercorbner, b_uppercorner) %>% group_by(mutant) %>% summarise(count = n()) 
# C2(1190) 60%, C3(323) 16%, C4(461) 24%, total(1974)

SPN_1.2_plot <- plot_tissue(sim, num_of_bins=100) +
  ggplot2::geom_rect(ggplot2::aes(
    xmin=b_lowercorbner[1], 
    xmax=b_uppercorner[1], 
    ymin=b_lowercorbner[2], 
    ymax=b_uppercorner[2]), 
    color= 'red', fill='white')

sim$sample_cells("SPN06_1.2", b_lowercorbner, b_uppercorner)

message("samples SPN06_1.1 and SPN06_1.2 extracted!")
sample_AB <- plot_tissue(sim)


#-------------------------------------------------------------------------------
#---------------------------- CHEMOTHERAPY ENTERS I ----------------------------
#-------------------------------------------------------------------------------
chemo1_start <- sim$get_clock()

sim$update_rates("C1", rates = c(growth = 0, death = 1))
sim$update_rates("C2", rates = c(growth = 0, death = 1))
sim$update_rates("C3", rates = c(growth = 0, death = 1))
sim$update_rates("C4", rates = c(growth = 0, death = 1))

sim$run_until(sim$var("C2") <= 100)

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
#--------------------------------- sample 2.1 ----------------------------------
bbox <- sim$search_sample(c("C5" = 2000), 50, 50)
c_lowercorbner <- bbox$lower_corner
c_uppercorner <- bbox$upper_corner

sim$get_cells(c_lowercorbner, c_uppercorner) %>% group_by(mutant) %>% summarise(count = n()) 
# C5(2216) 100%, total(2216)

SPN_2.1_plot <- plot_tissue(sim, num_of_bins=100) +
  ggplot2::geom_rect(ggplot2::aes(
    xmin=c_lowercorbner[1], 
    xmax=c_uppercorner[1], 
    ymin=c_lowercorbner[2], 
    ymax=c_uppercorner[2]), 
    color= 'red', fill='white')

sim$sample_cells("SPN06_2.1", c_lowercorbner, c_uppercorner)


#-------------------------------------------------------------------------------
#----------------------------- CHEMOTHERAPY ENTERS II --------------------------
#-------------------------------------------------------------------------------
chemo2_start <- sim$get_clock()

sim$update_rates("C5", rates = c(growth = 0, death = 1))

sim$run_until(sim$var("C5") <= 100)

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

#--------------------------------- sample 3.1 ----------------------------------
bbox <- sim$search_sample(c("C6" = 2000), 50, 50)
d_lowercorbner <- bbox$lower_corner
d_uppercorner <- bbox$upper_corner

sim$get_cells(d_lowercorbner, d_uppercorner) %>% group_by(mutant) %>% summarise(count = n()) 
# C6(2195) 100%, total(2195)

SPN_3.1_plot <- plot_tissue(sim, num_of_bins=100) +
  ggplot2::geom_rect(ggplot2::aes(
    xmin=d_lowercorbner[1], 
    xmax=d_uppercorner[1], 
    ymin=d_lowercorbner[2], 
    ymax=d_uppercorner[2]), 
    color= 'red', fill='white')

sim$sample_cells("SPN06_3.1", d_lowercorbner, d_uppercorner)


#--------------------------------- sample 3.2 ----------------------------------
e_lowercorbner <- as.integer(c(375, 390))
e_uppercorner <- as.integer(c(424, 439))

sim$get_cells(e_lowercorbner, e_uppercorner) %>% group_by(mutant) %>% summarise(count = n()) 
# C6(371) 18%, C7(1637) 82%, total(2008)

SPN_3.2_plot <- plot_tissue(sim, num_of_bins=100) +
  ggplot2::geom_rect(ggplot2::aes(
    xmin=e_lowercorbner[1], 
    xmax=e_uppercorner[1], 
    ymin=e_lowercorbner[2], 
    ymax=e_uppercorner[2]), 
    color= 'red', fill='white')

sim$sample_cells("SPN06_3.2", e_lowercorbner, e_uppercorner)

message("samples SPN06_3.1 and SPN06_3.2 extracted!")
#sample_DE <- plot_tissue(sim)

message("SIMUATION DONE!")


#-------------------------------------------------------------------------------
#----------------------------------- Save forest -------------------------------
#-------------------------------------------------------------------------------
forest <- sim$get_samples_forest()
forest$save(paste0(outdir, "sample_forest.sff"))

#-------------------------------------------------------------------------------
#----------------------------------- Save chemotherapy timing ------------------
#-------------------------------------------------------------------------------
chemo_timing <- list(
  chemo1_start = chemo1_start, # round(chemo1_start, 0)
  chemo1_end = chemo1_end,     # round(chemo1_end, 0)
  chemo2_start = chemo2_start, # round(chemo2_start, 0)
  chemo2_end = chemo2_end      # round(chemo2_end, 0)
)
saveRDS(chemo_timing, paste0(outdir, "chemo_timing.rds"))
message("DONE!, everything safe and sound")


#-------------------------------------------------------------------------------
#----------------------------------- PLOTS -------------------------------------
#-------------------------------------------------------------------------------
#plot_state(sim)
plot_tissue(sim)
#plot_muller(sim)




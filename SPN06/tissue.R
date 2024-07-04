
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

sim <- SpatialSimulation(name = "SPN06", width = 1000, height = 1000, save_snapshots = FALSE, seed = 3)
sim$history_delta <- 1

#-------------------------------------------------------------------------------
#---------------------- Rising of Clones I : [C1, C2, C3, C4] ------------------
#-------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++ clone 01 ++++++++++++++++++++++++++++++++++++
sim$add_mutant(name = "C1", growth_rates = 0.1, death_rates = 0.01) # add species
sim$place_cell("C1", 500, 500) # displace the initial cell in the tissue
sim$run_up_to_size("C1", 1000)


#+++++++++++++++++++++++++++++++++ clone 02 ++++++++++++++++++++++++++++++++++++
sim$update_rates("C1", rates = c(growth = 0.1, death = 0.1))
sim$add_mutant(name = "C2", growth_rates = 0.5, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C1"), "C2")
sim$run_up_to_size("C2", 1000)


#+++++++++++++++++++++++++++++++++ clone 03 ++++++++++++++++++++++++++++++++++++
sim$add_mutant(name = "C3", growth_rates = 1, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C3")
sim$run_up_to_size("C3",1000)


#+++++++++++++++++++++++++++++++++ clone 04 ++++++++++++++++++++++++++++++++++++
sim$add_mutant(name = "C4", growth_rates = 2, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("C2"), "C4")
sim$run_up_to_size("C4",3000)


#-------------------------------------------------------------------------------
#------------------------------ Sampling [A, B] --------------------------------
#-------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++ sample A ++++++++++++++++++++++++++++++++++++
bbox = sim$search_sample(c("C2" = 30, "C3" = 270), nw = 20, nh = 20)

'
plot_tissue(sim) + 
  geom_rect(
    xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1], 
    ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2], 
    fill = NA, color = "black"
  )
'

sim$sample_cells("A", bbox$lower_corner, bbox$upper_corner)


#+++++++++++++++++++++++++++++++++ sample B ++++++++++++++++++++++++++++++++++++
bbox = sim$search_sample(c("C2" = 45, "C3" = 15, "C4" = 15), nw = 30, nh = 30)

'
plot_tissue(sim) + 
  geom_rect(
    xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1], 
    ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2], 
    fill = NA, color = "black"
  )
'

sim$sample_cells("B", bbox$lower_corner, bbox$upper_corner)


p1 <- plot_tissue(sim)


#-------------------------------------------------------------------------------
#------------------------------ Treatment ENTERS I -----------------------------
#-------------------------------------------------------------------------------

chemo1_start <- sim$get_clock()

sim$update_rates("C1", rates = c(growth = 0, death = 3))
sim$update_rates("C2", rates = c(growth = 0, death = 3))
sim$update_rates("C3", rates = c(growth = 0, death = 3))
sim$update_rates("C4", rates = c(growth = 0, death = 3))

while ((sim$get_cells() %>% dplyr::filter(mutant == "C2") %>% nrow()) > 100) {
  sim$run_up_to_time(sim$get_clock() + 0.1)
}

chemo1_end <- sim$get_clock()


#-------------------------------------------------------------------------------
#--------------------------- Rising of Clones II : [C5] ------------------------
#-------------------------------------------------------------------------------

sim$add_mutant(name = "C5", growth_rates = 1, death_rates = 0.01)

sim$mutate_progeny(
  as.data.frame((sim$get_cells() %>% filter(mutant == 'C2') %>% arrange(position_x))[1,]), 
  "C5"
  )

sim$run_up_to_size("C5", 1e4)


#-------------------------------------------------------------------------------
#-------------------------------- Sampling [C] ---------------------------------
#-------------------------------------------------------------------------------

# sample C
bbox = sim$search_sample(c("C5" = 320), nw = 20, nh = 20)

'
plot_tissue(sim) + 
  geom_rect(
    xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1], 
    ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2], 
    fill = NA, color = "black"
  )
'

sim$sample_cells("C", bbox$lower_corner, bbox$upper_corner)

p2 <- plot_tissue(sim)


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
#------------------------- Rising of Clones III : [C6] -------------------------
#-------------------------------------------------------------------------------

sim$add_mutant(name = "C6", growth_rates = 1, death_rates = 0.01)

sim$mutate_progeny(
  as.data.frame((sim$get_cells() %>% filter(mutant == 'C5') %>% arrange(position_x))[1,]), 
  "C6"
)
sim$run_up_to_size("C6", 1e4)


#-------------------------------------------------------------------------------
#----------------------------- Sampling [D, E] ---------------------------------
#-------------------------------------------------------------------------------

bbox_width <- 20

#--------------------------------- sample D ++++++++++++++++++++++++++++++++++++
bbox1_p <- c(340, 470)
bbox1_q <- bbox1_p + bbox_width

#+++++++++++++++++++++++++++++++++ sample E ++++++++++++++++++++++++++++++++++++
bbox2_p <- c(420, 430)
bbox2_q <- bbox2_p + bbox_width

#+++++++++++++++++++++++++++++++++ View the boxes ++++++++++++++++++++++++++++++
'
plot_tissue(sim) +
  geom_rect(xmin = bbox1_p[1], xmax = bbox1_q[1],
            ymin = bbox1_p[2], ymax = bbox1_q[2],
            fill = NA, color = "black") +
  geom_rect(xmin = bbox2_p[1], xmax = bbox2_q[1],
            ymin = bbox2_p[2], ymax = bbox2_q[2],
            fill = NA, color = "black")
'

#+++++++++++++++++++++++++++++++++ sampling ++++++++++++++++++++++++++++++++++++
# sample D
sim$sample_cells("D", bottom_left = bbox1_p, top_right = bbox1_q)
# sample E
sim$sample_cells("E", bottom_left = bbox2_p, top_right = bbox2_q)


p3 <- plot_tissue(sim)
p4 <- plot_state(sim)
p5 <- plot_muller(sim)
p6 <- plot_timeseries(sim)

forest <- sim$get_samples_forest()

p7 <- plot_forest(forest) %>% annotate_forest(forest)


design <- "abcd
           eggg
           fggg"

p <- patchwork::wrap_plots(
  a = p1, 
  b = p2, 
  c = p3, 
  d = p4, 
  e = p5, 
  f = p6, 
  g = p7,
  guides = "auto", 
  design = design
)


#-------------------------------------------------------------------------------
#----------------------------------- Save --------------------------------------
#-------------------------------------------------------------------------------

#setwd("/u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06")

curr_dir <- getwd()

ggsave(filename = paste(curr_dir, "/plots/plot_01.png", sep = ""), plot = p, width = 16, height = 10, dpi = 300)

forest$save( paste(curr_dir, "/data/forest.sff", sep = "") )

chemo_timing <- list(
  chemo1_start = chemo1_start, # round(chemo1_start, 0)
  chemo1_end = chemo1_end,     # round(chemo1_end, 0)
  chemo2_start = chemo2_start, # round(chemo2_start, 0)
  chemo2_end = chemo2_end      # round(chemo2_end, 0)
)

saveRDS( chemo_timing, paste(curr_dir, "/data/chemo_timing.rds", sep = "") )









#devtools::install_github("caravagnalab/rRACES")
#setwd("/Users/azadsadr/Documents")
#install.packages("patchwork")

library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)

curr_dir <- getwd()

#readRDS( paste(curr_dir, "/data/chemo_timing.rds", sep = "") )

message("start loading forest...")
samples_forest <- load_samples_forest( paste(curr_dir, "/data/forest.sff", sep = "") )
message("finished loading forest!\n")

message("start loading phylo forest...")
phylo_forest <- load_phylogenetic_forest( paste(curr_dir, "/data/phylo_forest.sff", sep = "") )
message("finished loading phylo forest!")



annot_forest <- plot_forest(samples_forest) %>%
  annotate_forest(phylo_forest,
                  samples = T,
                  MRCAs = T,
                  exposures = T,
                  drivers=T,
                  add_driver_label = T)


p_exposure_timeline <- plot_exposure_timeline(phylo_forest)

labels <- get_relevant_branches(samples_forest)
sticks <- plot_sticks(samples_forest, labels)

message("so far good 1")

pl <- annot_forest + sticks + plot_layout(nrow = 3, design = 'A\nA\nB\nB\nC')

message("so far good 2")

ggsave( paste(curr_dir, "/plots/SPN06_mutations2222222.png", sep = "") , plot = pl, width = 210, height = 297, units = "mm", dpi = 300)

message("so far good 3")

ggsave( paste(curr_dir, "/plots/SPN06_exposure2222222.png", sep = "") , plot = p_exposure_timeline, width = 210, height = 297, units = "mm", dpi = 300)



#-------------------------------------------------------------------------------






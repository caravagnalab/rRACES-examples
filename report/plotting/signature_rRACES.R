#setwd('/orfeo/cephfs/scratch/cdslab/shared/races/data_for_report/SPN03/')

library(rRACES)
library(ggalluvial)
library(dplyr)
library(ggplot2)
library(patchwork)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/REPO_UPDATED/rRACES-examples/report/plotting/utils.R")




plot_exposure_evolution <- function(sample_forest,phylo_forest,snapshot){
  simulation <- rRACES::recover_simulation(snapshot)
  color_map_clones <- get_clone_map(sample_forest)
  muller <- plot_muller(simulation,color_map = color_map_clones)
  
  signature <- phylo_forest$get_exposures()
  df_sign <- get_exposure_ends(phylo_forest)
  start_time <- 0
  end_time <- max(df_sign$end_time)
  
  df_final <- tibble()
  for (i in seq(start_time, end_time, 2)){
    for (irow in 1:nrow(df_sign)){
      row <- df_sign[irow,]
      if (i %in% seq(row$time, row$end_time)){
        tmp <- tibble(time = i, signature = row$signature, exposure = row$exposure, type = row$type)
        df_final <- bind_rows(df_final, tmp)
      }
    }
  }
  #
  #
  sign_1 <- ggplot(df_final, aes(x = factor(time), stratum = signature, alluvium = signature, y = exposure)) +
    geom_alluvium(aes(fill = signature), width = 0.001, curve_type = 'linear', decreasing = NA, knot.pos = 0) +
    scale_x_discrete(breaks = seq(0, 300, 20)) +
    labs(legend = 'Signature',
         x = "",
         y = "Exposure") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt")) + facet_grid(type~.)
  
  
  sign_muller <- sign_1 + muller + plot_layout(nrow = 2)
  
  
  
  
  sign_tree  <- plot_forest(sample_forest,color_map = color_map_clones) +
    df_final %>%
    group_by(time) %>%
    reframe(exposure = exposure/sum(exposure), across(everything())) %>%
    ggplot()+
    geom_bar(aes(x = time, y = exposure, fill = signature), stat = "identity", width = 10) +
    ylab('Exposure') +
    xlab('') +
    scale_x_reverse() +
    #scale_fill_manual('Signature', values = signatures_palette(forest_muts,55))+
    theme_minimal() +
    coord_flip() + theme(axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         plot.margin = margin(0, 0, 0, 0, "pt")) + facet_grid(.~type)
  return(list("sign_tree"=sign_tree,"sign_muller"=sign_muller))
}


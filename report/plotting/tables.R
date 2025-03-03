library(rRACES)
library(tidyverse)

samples_table <- function(snapshot, forest) {
    sim <- rRACES::recover_simulation(snapshot)
    sample_forest <- load_samples_forest(forest)
    info = sim$get_samples_info() ## requested from either the simulation recovery or as saved table

    nodes = sample_forest$get_nodes()
    clones = nodes %>% 
      dplyr::filter(!is.na(sample)) %>% 
      dplyr::group_by(sample, mutant) %>% 
      dplyr::pull(mutant) %>% 
      unique()
    clones_of_origin = nodes %>%
      dplyr::filter(!is.na(sample)) %>% 
      dplyr::group_by(sample, mutant) %>% 
      # dplyr::mutate(mutant = gsub(" ", "_", mutant)) %>% 
      dplyr::count(mutant) %>% 
      tidyr::pivot_wider(values_from = n, names_from = mutant, values_fill = 0) %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(Sample_Type = sum(c_across(clones) == 0)) %>% 
      dplyr::mutate(Sample_Type = ifelse(Sample_Type == (length(clones)-1), "Monoclonal", "Polyclonal")) %>% 
      dplyr::mutate(Total_Cells = sum(c_across(clones), na.rm = T)) %>% 
      dplyr::mutate(across(all_of(clones), ~ round(.x/Total_Cells,2), .names = "{.col} proportion"))
    
    samples_tb = dplyr::full_join(info, clones_of_origin, by = join_by("name" == "sample")) %>% 
      select(!c("xmin","xmax","ymin","ymax","id","tumour_cells","tumour_cells_in_bbox")) %>% 
      rename(Sample_ID=name) %>% 
      rename(Samping_Time=time) %>% 
      mutate(Samping_Time=round(Samping_Time,2)) %>% 
      arrange(Sample_ID)
      


    return(samples_tb)
}



library(ProCESS)
library(tidyverse)

# Sys.setenv(DOWNLOAD_STATIC_LIBV8=1)
# install.packages('V8')
install.packages('gt')
library(gridtext)
library(gt)

ratings_in_time = function(sim, time_points) {
  options(pillar.sigfig = 5, digits = 5)
  history_rates = sim$get_rates_update_history() %>% 
        filter(time %in% time_points) %>% 
        tidyr::pivot_wider(values_from = rate, names_from = event) %>%
        dplyr::rename(death_rate = death, growth_rate = growth) %>%
        dplyr::relocate(death_rate, .after = growth_rate) %>% 
        dplyr::arrange(time)

    history_counts = sim$get_count_history() %>% 
        dplyr::filter(time %in% time_points) %>%
        dplyr::rename(n_cells = count)

    tb = dplyr::right_join(history_rates, history_counts, 
        by = join_by("time" == "time", "mutant" == "mutant", "epistate" == "epistate")) %>%
        dplyr::filter(!is.na(growth_rate))
    # adding the final number of cells and growth rates of present clones 
    
    rn = sim$get_clock()
    rn_counts = sim$get_counts()
    rn_clones = rn_counts %>% 
        dplyr::filter(counts != 0) %>%
        dplyr::pull(mutant)
    rn_rates = lapply(rn_clones, function(x) {
        dplyr::tibble(mutant = x, 
            growth_rate = sim$get_rates(x)$growth, 
            death_rate = sim$get_rates(x)$death)
    }) %>% dplyr::bind_rows()
    rn_info = dplyr::right_join(rn_counts, rn_rates, by = "mutant") %>%
        dplyr::mutate(time = rn) %>%
        dplyr::relocate(time, .before = mutant) %>%
        dplyr::relocate(counts, .after = death_rate) %>%
        dplyr::mutate(overall = NULL) %>%
        dplyr::rename(n_cells = counts)

    tb = dplyr::bind_rows(tb, rn_info)
    
    moment_zero = history_rates[1, ] %>%
        dplyr::mutate(n_cells = 0)
    
    tb = dplyr::bind_rows(moment_zero, tb)
    tb = tb %>% 
      dplyr::group_by(time) %>% 
      dplyr::group_split()
    
    tb = lapply(1:length(tb), function(x) {
      tt = as.numeric(x) - 1
      tb[[x]] %>% 
        dplyr::mutate(t = case_when(as.numeric(time) == rn ~ paste0("t~", tt, "~", " (last time point)"), .default = paste0("t~", tt, "~")))
    }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::relocate(t, .before = time)
    
    tt = tb %>%
      dplyr::group_by(t) %>% 
      gt::gt(row_group_as_column = T, process_md = T) %>%
      # gt::gt(groupname_col = "t", process_md = T) %>%
      gt::tab_stubhead(label = md("**Time point**")) %>% 
      gt::tab_header(title = md("**Rates and number of cells of each clone during tissue growth**"), 
                     subtitle = md("At each time point the growth and/or death rate of one clone changed")) %>% 
      gt::cols_label(mutant = md("Mutant"), 
                     epistate = md("Epistate"), 
                     growth_rate = md("Growth rate"), 
                     death_rate = md("Death rate"), 
                     n_cells = md("number of cells"), 
                     time = md("Time")) %>% 
      gt::tab_footnote(locations = cells_column_labels(time), footnote = "Simulation time")
    tt = gt::as_gtable(tt, text_grob = gridtext::richtext_grob, plot = T)
    return(tt)
}


samples_table <- function(sim, forest) {
    info = sim$get_samples_info()
    
    nodes = forest$get_nodes()
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
      dplyr::mutate(sample_type = sum(c_across(clones) == 0)) %>% 
      dplyr::mutate(sample_type = ifelse(sample_type == (length(clones)-1), "Monoclonal", "Polyclonal"))
    
      # dplyr::mutate(total_cells = sum(c_across(clones), na.rm = T)) 
    
    info = dplyr::full_join(info, clones_of_origin, by = join_by("name" == "sample"))

    info = info %>% 
      dplyr::group_by(time) %>% 
      dplyr::group_split() 
    info = lapply(1:length(info), function(x) {
        oo = info[[x]] %>% 
          dplyr::mutate(t = paste0("t~", x, "~"))
      }) %>% 
      bind_rows()
    
    samples_tb = info %>% 
      dplyr::select(t, name, tumour_cells, time, all_of(clones), sample_type) %>% 
      dplyr::group_by(t) %>% 
      dplyr::arrange(name) %>% 
      gt::gt(row_group_as_column = T, process_md = T) %>% 
      gt::tab_stubhead(label = md("**Time point**")) %>% 
      gt::tab_header(title = md("**Tumour samples**")) %>% 
      gt::cols_label(name = md("Sample"), 
                     tumour_cells = md("Number of sampled cells"), 
                     time = md("Time of sampling"), 
                     sample_type = md("Sample type")) %>% 
      gt::tab_footnote(locations = cells_stubhead(), footnote = "Sampling time") %>% 
      gt::tab_footnote(locations = cells_column_labels(time), footnote = "Simulation time") %>% 
      gt::data_color(columns = sample_type, direction = "column", apply_to = "text", palette = c("Polyclonal" = "#C40C0C", "Clonal" = "#333A73"))
    
    samples_tb = gt::as_gtable(samples_tb, text_grob = gridtext::richtext_grob, plot = T)
    return(samples_tb)
}



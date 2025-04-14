my_muller_plot = function (simulation) 
{
  stopifnot(inherits(simulation, "Rcpp_Simulation"))
  df_populations <- simulation$get_count_history() %>% dplyr::as_tibble() %>% 
    dplyr::mutate(Identity = paste0(.data$mutant, .data$epistate)) %>% 
    dplyr::rename(Generation = .data$time, Population = .data$count) %>% 
    dplyr::select(.data$Generation, .data$Identity, .data$Population)
  df_populations = df_populations %>% mutate(Population = log(Population, base=10)) %>%
    mutate(Population = ifelse(Population == -Inf, 0, Population))
  df_edges <- simulation$get_lineage_graph() %>% dplyr::distinct(.data$ancestor, 
                                                                 .data$progeny) %>% dplyr::rename(Parent = .data$ancestor, 
                                                                                                  Identity = .data$progeny) %>% dplyr::select(.data$Parent, 
                                                                                                                                              .data$Identity)
  df_edges <- ProCESS:::collapse_loops(df_edges)
  max_tumour_size <- df_populations %>% dplyr::group_by(.data$Generation) %>% 
    dplyr::summarise(Population = sum(.data$Population)) %>% 
    dplyr::pull(.data$Population) %>% max()
  #max_tumour_size = log(max_tumour_size, base=10)
  max_tumour_size <- max_tumour_size * 1.05
  wt_dynamics <- df_populations %>% dplyr::group_by(.data$Generation) %>% 
    dplyr::summarise(Population = sum(.data$Population)) %>% 
    dplyr::mutate(Identity = "Wild-type", Population = max_tumour_size - 
                    .data$Population)
  t_wt_dynamics <- dplyr::bind_rows(wt_dynamics, df_populations)
  color_map <- ProCESS:::get_species_colors(simulation$get_species())
  suppressWarnings({
    muller_df <- ggmuller::get_Muller_df(df_edges, t_wt_dynamics)
    plot <- ggmuller::Muller_pop_plot(muller_df, add_legend = TRUE, 
                                      palette = c(`Wild-type` = "gainsboro", color_map)) + 
      ProCESS:::my_theme() + ggplot2::guides(fill = ggplot2::guide_legend("Species"))
  })
  plot
}
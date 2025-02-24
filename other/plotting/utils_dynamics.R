get_timepoint_info <- function(sim, t){
  tissue <- plot_tissue(sim)+ggtitle(t)
  state <- plot_state(sim)+ggtitle(t)
  muller <- plot_muller(sim)+ggtitle(t)
  fire <- plot_firings(sim)+ggtitle(t)
  timeserie <- plot_timeseries(sim)+ggtitle(t)
  time <- sim$get_clock()
  label <- t
  return(list(tissue=tissue,
              state=state,
              muller=muller,
              fire=fire,
              timeserie = timeserie,
              time=time,
              label = label))
}

get_sampling_info <- function(sim, sample_name, bbox){
  box_data <- tibble(
    sample = unname(sample_name), 
    x_min = bbox$lower_corner[1],
    x_max = bbox$upper_corner[1],
    y_min = bbox$lower_corner[2],
    y_max = bbox$upper_corner[2])
  
  plot <- plot_tissue(sim) + 
    ggtitle(sample_name) +
    geom_rect(data = box_data, aes(
      xmin = x_min,
      xmax = x_max,
      ymin = y_min,
      ymax = y_max),
      color = 'indianred3',
      fill = NA, 
      inherit.aes = FALSE) +
    geom_text_repel(data = box_data,
                    aes(y = y_max, x = x_max, label = sample), size=3, inherit.aes = FALSE, show.legend = F, 
                    force = 2, min.segment.length = 0,
                    position = position_nudge_repel(x = -0.1, y = 0.05)) 
  time <-  sim$get_clock()
  label <- sample_name
  return(list(plot = plot,
              time = time, 
              label = label))
}

get_clone_info <- function(sim, forest){
  info <- sim$get_samples_info()
  nodes <- forest$get_nodes()
  clones <- nodes %>% 
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    dplyr::pull(mutant) %>% 
    unique()
  clones_of_origin <- nodes %>%
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    dplyr::count(mutant) %>% 
    tidyr::pivot_wider(values_from = n, names_from = mutant, values_fill = 0) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(sample_type = sum(c_across(clones) == 0)) %>% 
    dplyr::mutate(sample_type = ifelse(sample_type == (length(clones)-1), "Monoclonal", "Polyclonal"))
  
  final_clones <- sim$get_counts() %>% filter(counts != 0) %>% pull(mutant)
  clones_of_origin <- clones_of_origin %>% tidyr::pivot_longer(cols = final_clones)
  
  return(clones_of_origin)
}

plot_clonal_prop <- function(info, sim){
  color_map <- rRACES:::get_species_colors(sim$get_species())
  pi_sampling <- lapply(unique(info$sample), FUN = function(s){
    t <- info %>% filter(sample == s) 
    ggplot2::ggplot(t) +
      ggplot2::geom_bar(stat = "identity",
                        ggplot2::aes(x = "",
                                     y = .data$value,
                                     fill = .data$name)) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::theme_void(base_size = 10) +
      ggplot2::scale_fill_manual(values = color_map) +
      ggplot2::theme(legend.position = "bottom")  +
      ggplot2::labs(fill = "Species") + 
      ggtitle(s)
  }) 
  
  return(pi_sampling)
}

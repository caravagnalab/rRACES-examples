library(ProCESS)
library(ggplot2)
library(patchwork)
library(ggrepel)

plot_tumour_dynamics <- function(snapshot,sample_forest){

  simulation <- ProCESS::recover_simulation(snapshot)
  simulation_info <- simulation$get_samples_info()
  color_map_clones <- get_clone_map(sample_forest)
  muller <- plot_muller(simulation,color_map = color_map_clones)
  timeseries <- plot_timeseries(simulation, color_map = color_map_clones)


  timing <- simulation_info %>% mutate(what = 'simulation time')
  sim_end <- max(timing$time)
  time_plot <- ggplot(timing, aes(x = time, y = what, label = name)) +
    geom_line() +
    geom_point() +
    ggrepel::geom_text_repel(direction = "y",
                             point.padding = 0.5,
                             hjust = 0,
                             box.padding = 1,
                             seed = 123)  +
    xlim(0, sim_end) +
    scale_y_discrete(name = "") +
    theme_minimal()


  samples_name <- simulation_info$name
  plot_sample <- lapply(samples_name, FUN = function(s){
    plot_tissue(simulation = simulation, before_sample = s,color_map = color_map_clones)  + ggtitle(s)
  })
  fixed_n_samples_per_row<-3
  n_row <- ceiling(length(samples_name)/fixed_n_samples_per_row)
  plot_sampling <- patchwork::wrap_plots(plot_sample, nrow = n_row, ncol = 3, guides = 'collect') & theme(legend.position = 'bottom')


  plot_dynamics <- time_plot + muller + timeseries + plot_layout(nrow = 3)
  return(list("plot_sampling" = plot_sampling, "plot_dynamics" = plot_dynamics))
}

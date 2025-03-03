library(rRACES)
library(ggplot2)
library(patchwork)
library(ggrepel)

plot_tumour_dynamics <- function(snapshot){

  simulation <- rRACES::recover_simulation(snapshot)
  simulation_info <- simulation$get_samples_info()

  muller <- plot_muller(simulation)
  timeseries <- plot_timeseries(simulation)


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
    plot_tissue(simulation = simulation, before_sample = s)  + ggtitle(s)
  })

  plot_sampling <- patchwork::wrap_plots(plot_sample, nrow = 1, guides = 'collect') & theme(legend.position = 'bottom')


  final_plot <- time_plot + muller + timeseries + plot_layout(nrow = 3)
  return(list("plot_sampling" = plot_sampling, "plot_dynamics" = final_plot))
}
snapshot<-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/races/SPN01"
plot <- plot_tumour_dynamics(snapshot)
saveRDS(file = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/races/plot_dynamics.rds",
        object = plot)
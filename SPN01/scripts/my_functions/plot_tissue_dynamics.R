plot_tissue_dynamics<-function(tissue_plots,state_plots,SPN_id){
  if (length(tissue_plots)!=length(state_plots)){
    stop("Tissue plots and timeseries plots must have the same lenght.")
  }
  layout_dinamycs = "
    AABBCCDD
    EEFFGGHH
    "
  layout_pairs <- c(
    area(t = 1, l = 1, b = 5, r = 5),
    area(t = 1, l = 4, b = 3, r = 5)
  )
  list_pairs_plots <- lapply(seq_along(rep(1,8)), function(x){
				           ggplot2::ggplot()
					     })
  x =  patchwork::wrap_plots(list_pairs_plots, design = layout_dinamycs)
  for (i in seq_along(tissue_plots)){
    t <- tissue_plots[[i]] + ggplot2::theme(legend.position = "none") 
    p <- state_plots[[i]]+ggplot2::theme(legend.position = "none") 
    list_pairs_plots[[i]] <- t + p +
      plot_layout(design = layout_pairs)
  }
  last <- length(state_plots)
  leg <- ggpubr::get_legend(state_plots[[last]])
  plot_legend <- ggpubr::as_ggplot(leg)
  x =  patchwork::wrap_plots(list_pairs_plots, design = layout_dinamycs) 
  final_plot <- x+ plot_annotation(title = "Tissue dinamycs")
  return(final_plot)
}





#' plotting_sample
#'
#' @param sim rRACES simulation object
#' @param samples_timing a list of vector. names are the time point of sampling, and each element is a vector with names (indicating if the sample is clonal or not) 
#' with name of the samples present at that time point
#' @param boxes a list of Rcpp_TissueRectangle objects from rRACES. Must have the names of the correspondent sample. 
#'
#' @return plot
#' @export
#'
#' @examples 
#' bbox = new(TissueRectangle, c(430, 350), 21, 21)
#' bbox2 = new(TissueRectangle, c(300, 510), 21, 21)
#' s = c("Clonal" = "A", "Polyclonal" = "B")
#' timing = list("T1" = s)
#' box = list("A" = bbox, "B" = bbox2)
#' plotting_sample(sim, samples_timing = timing, boxes = box)

library(cli)
library(english)
library(ggrepel)

plotting_sample = function(sim, samples_timing, boxes) {
  
  color_type = c("Polyclonal" = "#C40C0C", "Clonal" = "#333A73")
  
  layout = "
    AABBCCDD
    "
  
  p_list <- lapply(english::ordinal(seq(1:4)), function(x) {ggplot() + ggtitle(paste(x, "sampling"))})
  
  for (t in seq_along(samples_timing)) {
    samples = samples_timing[[t]]
    print(samples)
    p_list[[t]] <- rRACES::plot_tissue(sim) 

    for (i in seq_along(samples)) {
      
      sample = samples[i]
      bbox = boxes[[sample]]
      
      if(class(bbox)[1] != "Rcpp_TissueRectangle") {
        cli::cli_abort("Must provide a list of {.cls Rcpp_TissueRectangle} objects as {.field boxes} argument")
      }
      
      # print(sample)
      sample_data = tibble(type = names(sample), 
                           sample = unname(sample), 
                           x_min = bbox$lower_corner[1],
                           x_max = bbox$upper_corner[1],
                           y_min = bbox$lower_corner[2],
                           y_max = bbox$upper_corner[2])
      
      p_list[[t]] =  p_list[[t]] +
        geom_rect(data = sample_data, aes(color = type, 
                                          xmin = x_min,
                                          xmax = x_max,
                                          ymin = y_min,
                                          ymax = y_max),
                  fill = NA, inherit.aes = FALSE) +
        scale_color_manual(values = color_type) +
        # geom_text_repel(data = sample_data, aes(label = sample, color = type), size=3, inherit.aes = FALSE, show.legend = F, check_overlap = T) +
        geom_text_repel(data = sample_data, aes(y = y_max, x = x_max, label = sample, color = type), size=3, inherit.aes = FALSE, show.legend = F, 
                        force = 2, min.segment.length = 0,
                         position = position_nudge_repel(x = -0.1, y = 0.05)) +
        ggtitle(paste(english::ordinal(t), "sampling")) +
        guides(color = guide_legend(nrow = 1, title = "Type of sample", title.position = "top")) +
        theme(legend.position = "bottom")
    }
  }
  
  p = wrap_plots(p_list, design = layout)+ plot_layout(guides = 'collect') & ggplot2::theme(legend.position = "none")
  #/ guide_area()  + plot_layout(guides = 'collect')
  return(p)
}



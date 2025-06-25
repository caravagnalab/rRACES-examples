
rm(list = ls())
options(bitmapType='cairo')
require(tidyverse)
library(optparse)

source("../../getters/process_getters.R")
source("utils/plot_utils.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN03'),
	            make_option(c("--coverages"), type = "character", default = '50, 100'),
		    make_option(c("--purities"), default = '0.3, 0.6, 0.9')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'
spn_id = opt$spn_id

cleaned <- gsub('^"|"$', '', opt$purities)
purity_list <- strsplit(cleaned, ",")[[1]]
PURITY <- trimws(purity_list)
cleaned <- gsub('^"|"$', '', opt$coverages)
coverage_list <- strsplit(cleaned, ",")[[1]]
COVERAGES <- trimws(coverage_list)
print(COVERAGES)
print(PURITY)

MUT_TYPES = c("INDEL", "SNV")

params_grid = expand.grid(COVERAGES, PURITY, MUT_TYPES)
colnames(params_grid) = c("coverage", "purity", "mut")

input_dir <-  paste0(data_dir,spn_id,"/validation/somatic/")

i = 1
df = lapply(1:nrow(params_grid), function(i) {
  coverage = params_grid[i,]$coverage
  purity = params_grid[i,]$purity
  mut_type = params_grid[i,]$mut
  combination = paste0(coverage, "x_", purity, "p")
  
  results_folder_path = file.path(input_dir, spn_id, combination, "allCaller", mut_type)
  file_name = file.path(results_folder_path, "metrics.rds")
  if (file.exists(file_name)) {
    dplyr::bind_cols(readRDS(file_name), params_grid[i,])    
  }
}) %>% do.call("bind_rows", .)


plot_rep = function(df) {
  # 1. Precision comparison with mut distinction via linetype
  p1 <- ggplot(df, aes(x = VAF_bin, y = precision, color = caller, linetype = mut)) +
    geom_line(aes(group = interaction(caller, mut)), size = 1) +
    geom_point(size = 2) +
    facet_grid(coverage ~ purity, labeller = label_both) +
    labs(title = "Precision by VAF Bin Across Coverage and Purity",
         x = "VAF Bin", y = "Precision",
         color = "Caller", linetype = "Mutation Type") +
    theme_bw() +
    scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
    scale_color_manual(values = method_colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 8)) +
    scale_y_continuous(limits = c(0, 1))
  
  # 2. Sensitivity/Recall comparison with mut distinction via linetype
  p2 <- ggplot(df, aes(x = VAF_bin, y = sensitivity, color = caller, linetype = mut)) +
    geom_line(aes(group = interaction(caller, mut)), size = 1) +
    geom_point(size = 2) +
    facet_grid(coverage ~ purity, labeller = label_both) +
    labs(title = "Sensitivity (Recall) by VAF Bin Across Coverage and Purity",
         x = "VAF Bin", y = "Sensitivity",
         color = "Caller", linetype = "Mutation Type") +
    theme_bw() +
    scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
    scale_color_manual(values = method_colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 8)) +
    scale_y_continuous(limits = c(0, 1))
  
  # 5. Combined metric plot with mut distinction via shape
  p5 <- ggplot(df, aes(x = sensitivity, y = precision, color = caller, shape = mut)) +
    geom_point(aes(size = vaf_correlation), alpha = 0.7) +
    facet_grid(coverage ~ purity, labeller = label_both) +
    labs(title = "Precision vs Sensitivity",
         x = "Sensitivity", y = "Precision",
         color = "Caller", size = "VAF Corr", shape = "Mutation Type") +
    scale_color_manual(values = method_colors) +
    theme_bw() +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5)
  
  # 6. Performance summary with mut as additional faceting dimension
  # Reshape data for multiple metrics
  df_long <- df %>%
    select(caller, coverage, purity, mut, VAF_bin, precision, sensitivity, 
           vaf_correlation) %>%
    tidyr::pivot_longer(cols = c(precision, sensitivity, vaf_correlation),
                        names_to = "metric", values_to = "value")
  
  # Alternative version of p6 with mut as fill instead of facet (if too many facets)
  p6 <- ggplot(df_long, aes(x = caller, y = value, fill = mut)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    facet_grid(metric ~ paste("Cov:", coverage, "Pur:", purity), 
               scales = "free_y", labeller = label_value) +
    labs(title = "Performance Metrics Distribution by Caller and Mutation Type",
         x = "Caller", y = "Metric Value", fill = "Mutation Type") +
    theme_bw() +
    scale_fill_manual(values = c("INDEL" = "steelblue", "SNV" = "darkorange")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 7))
  
  design <- "
  AB
  CD
  "
  
  # Generate descriptive title and subtitle for the report
  title <- paste0(spn_id, " - somatic mutations")
  
  # Combine all plots into a comprehensive report using patchwork
  # free() function allows each plot to maintain its own scales
  report_plot <- patchwork::free(p1) + patchwork::free(p2) +
    patchwork::free(p5) + patchwork::free(p6) +
    patchwork::plot_layout(design = design) +
    patchwork::plot_annotation(title) & 
    ggplot2::theme(text = ggplot2::element_text(size = 12))
  
  report_plot
}

final_report = plot_rep(df)
ggsave(paste0(input_dir,"example_report.pdf"), plot = final_report, width = 18, height = 14)


rm(list = ls())
options(bitmapType='cairo')
require(tidyverse)
library(optparse)

source("../../getters/process_getters.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN03'))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'
spn_id = opt$spn_id

spn_id = "SPN04"

COVERAGES = c("50", "100", "150", "200")
PURITY = c("0.3", "0.6", "0.9")
CALLERS = c("mutect2", "strelka", "freebayes")
MUT_TYPES = c("INDEL", "SNV")
SAMPLES = get_sample_names(spn = spn_id)

params_grid = expand.grid(COVERAGES, PURITY, CALLERS, MUT_TYPES, SAMPLES)
colnames(params_grid) = c("coverage", "purity", "caller", "mut", "sample")

input_dir <-  paste0(data_dir,spn_id,"/validation/somatic/")
outdir <- paste0(data_dir,spn_id,"/validation/somatic/report")
dir.create(outdir, recursive = T, showWarnings = F)

metrics_df = dplyr::tibble()
vaf_df = dplyr::tibble()

print_missing = FALSE
print_avail = TRUE

for (i in 1:nrow(params_grid)) {
  
  params = params_grid[i,]
  
  coverage = params$coverage
  purity = params$purity
  caller = params$caller
  mut_type = params$mut
  sample_id = params$sample
  combination = paste0(coverage, "x_", purity, "p")
  
  folder_path = file.path(input_dir, spn_id, combination, caller, sample_id, mut_type)
  
  all_metrics_path = file.path(folder_path, "metrics.rds")
  
  if (file.exists(all_metrics_path)) {
    all_metrics = readRDS(all_metrics_path)
    if ("chr_caller" %in% colnames(all_metrics$vaf_comparison)) {
      if (print_avail) print(paste0(all_metrics_path, " avail"))
      
      metrics_df = dplyr::bind_rows(metrics_df, cbind(all_metrics$report_metrics, params))
      vaf_df = dplyr::bind_rows(vaf_df, cbind(all_metrics$vaf_comparison, params))  
    } else {
      if (print_missing) print(paste0(all_metrics_path, " not avail"))
    }
  } else {
    if (print_missing) print(paste0(all_metrics_path, " not avail"))
  }
}


metrics_df %>% 
  dplyr::filter(mut == "SNV") %>% 
  ggplot(mapping = aes(x = caller, y = value, fill = name)) +
  geom_col(position = "dodge") +
  facet_grid(purity~coverage) +
  theme_bw()

metrics_df %>% 
  dplyr::mutate(purity = as.numeric(as.character(purity))) %>% 
  dplyr::filter(mut == "SNV", Mutations == "Only Pass") %>% 
  ggplot(mapping = aes(x = coverage, y = value, col = caller, alpha = purity)) +
  geom_point() +
  facet_grid(sample~name) +
  theme_bw()

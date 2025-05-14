
rm(list = ls())
require(tidyverse)
source("utils/utils.R")
source("utils/plot_utils.R")

read_all_metrics = function(spn, mut_type, sample_id) {
  all_metrics = dplyr::tibble()
  all_vaf_comparison = dplyr::tibble()
  
  if (!dir.exists(spn)) {stop("No directory exist for the given spn")}
  purity_cov_combinations = list.files(spn)
  if (length(purity_cov_combinations) == 0) {stop("No results are parsed for the given SPN")}
  for (c in purity_cov_combinations) {
    message(paste0("Parsing ", c, " combination"))
    purity = as.numeric(str_split(c, "_")[[1]][2])
    coverage = as.numeric(str_replace(str_split(c, "_")[[1]][4], "x", ""))
    callers = list.files(file.path(spn, c))
    callers = callers[callers!="races"]
    for (call in callers) {
      metrics_path = file.path(spn,c,call,sample_id,mut_type,"metrics.rds")
      if (file.exists(metrics_path)) {
        message("    Got ", call, " res")
        metrics = readRDS(metrics_path)
        
        all_metrics = dplyr::bind_rows(
          all_metrics, 
          metrics$report_metrics %>% 
            dplyr::mutate(spn=spn, caller=call, mut_type=mut_type, purity=purity, coverage=coverage, sample_id=sample_id)
        )
        
        all_vaf_comparison = dplyr::bind_rows(
          all_vaf_comparison, 
          metrics$vaf_comparison %>% 
            dplyr::mutate(spn=spn, caller=call, mut_type=mut_type, purity=purity, coverage=coverage, sample_id=sample_id)
        )
        
      } else {
        message("    No results for ", call)
      }
    }
  }
  
  list(vaf_comparison=all_vaf_comparison, metrics=all_metrics)
}

# sample_id = "SPN01_1.2"
# mut_type = "INDEL"
# spn = "SPN01"

df_vaf_errors = dplyr::tibble()
df_all_metrics = dplyr::tibble()
for (sample_id in c("SPN01_1.1", "SPN01_1.2", "SPN01_1.3")) {
  metrics_snv = read_all_metrics("SPN01", "SNV", sample_id)
  metrics_indel = read_all_metrics("SPN01", "INDEL", sample_id)
  df_vaf_errors = dplyr::bind_rows(df_vaf_errors, metrics_snv$vaf_comparison, metrics_indel$vaf_comparison)
  df_all_metrics = dplyr::bind_rows(df_all_metrics, metrics_snv$metrics, metrics_indel$metrics)
}

# VAF error plot
alpha = .05
df_vaf_errors %>% 
  dplyr::mutate(coverage = factor(coverage)) %>% 
  dplyr::mutate(VAF_error = 100 * abs(VAF_difference) / VAF_truth) %>% 
  dplyr::group_by(caller, sample_id, coverage, purity, mut_type) %>% 
  dplyr::summarise(mean_err = mean(VAF_error), low_err = stats::quantile(VAF_error, alpha), high_err = stats::quantile(VAF_error, 1-alpha)) %>% 
  ggplot(mapping = aes(x = purity, y=mean_err, ymin=low_err, ymax=high_err, col=caller, linetype=coverage)) +
  geom_pointrange(position = position_dodge(width = .015)) +
  geom_line(position = position_dodge(width = .015)) +
  theme_bw() +
  labs(x = "Purity", y="VAF error (%)", col = "Caller") +
  facet_grid(mut_type~sample_id)


my_colors = c(
  "freebayes SNV" = "#a6cee3",
  "freebayes INDEL" = "#1f78b4",
  "strelka SNV" = "#b2df8a",
  "strelka INDEL" = "#33a02c",
  "mutect2 SNV" = "#fb9a99",
  "mutect2 INDEL" = "#e31a1c"
)

df_vaf_errors %>% 
  dplyr::mutate(caller = paste0(caller," ", mut_type)) %>% 
  dplyr::mutate(coverage = factor(coverage)) %>% 
  dplyr::mutate(VAF_error = 100 * abs(VAF_difference) / VAF_truth) %>% 
  dplyr::group_by(caller, sample_id, coverage, purity) %>% 
  dplyr::summarise(mean_err = mean(VAF_error), low_err = stats::quantile(VAF_error, alpha), high_err = stats::quantile(VAF_error, 1-alpha)) %>% 
  ggplot(mapping = aes(x = purity, y=mean_err, ymin=low_err, ymax=high_err, col=caller, linetype=coverage)) +
  geom_pointrange(position = position_dodge(width = .015)) +
  geom_line(position = position_dodge(width = .015)) +
  theme_bw() +
  labs(x = "Purity", y="VAF error (%)", col = "Caller") +
  facet_grid(coverage~sample_id, scales = "free_y") +
  scale_color_manual(values = my_colors)

df_vaf_errors %>% 
  dplyr::mutate(caller = paste0(caller," ", mut_type)) %>% 
  dplyr::mutate(coverage = factor(coverage)) %>% 
  dplyr::mutate(VAF_error = 100 * abs(VAF_difference) / VAF_truth) %>% 
  ggplot(mapping = aes(x = factor(purity), y=VAF_error, group=paste0(purity, caller), fill=caller)) +
  geom_boxplot(outliers = F) +
  theme_bw() +
  labs(x = "Purity", y="VAF error (%)", col = "Caller") +
  facet_grid(sample_id~coverage, scales = "free_y") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
  scale_y_continuous(trans="log10")


# Metrics plot
df_all_metrics %>% 
  dplyr::mutate(coverage = factor(coverage)) %>% 
  dplyr::mutate(caller = paste0(caller," ", mut_type)) %>% 
  dplyr::filter(Mutations != "All") %>% 
  ggplot(mapping = aes(x=purity, y=value, col=caller, linetype=coverage)) +
  geom_point() +
  geom_line() +
  facet_grid(sample_id~name) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_bw()

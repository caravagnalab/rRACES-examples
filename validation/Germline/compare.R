library(tidyverse)
library(vcfR)
library(optparse)
library(caret)
library(dplyr)

setwd('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/ProCESS-examples/validation/Germline/')
source('vcf_parser.R')
source('utils.R')

option_list <- list( 
  make_option(c("-i", "--input"), type="character", default='/orfeo/LTS/LADE/LT_storage/lvaleriani/races/validation_data', help="path to input data"),
  make_option(c("-s", "--SPN"), type="character", default='SPN01', help="SPN name"),
  make_option(c("-t", "--tool"), type="character", default='freebayes', help="variant calling tool"),
  make_option(c("-o", "--output"), type="character", default='/orfeo/LTS/LADE/LT_storage/lvaleriani/races/validation_data', help="path to output directory")
)
to_parse = TRUE

param <- parse_args(OptionParser(option_list=option_list))
dir <- param$input
spn <- param$SPN
tool <- param$tool
output <- param$output

vcf_list = lapply(1:22, FUN = function(chr){
  vcf_file = paste0(dir,'/', spn, '/normal_sample/', tool, '/vcf/chr', chr, '_normal_sample.', tool, '.vcf.gz')
  rds_file = paste0(dir,'/', spn, '/normal_sample/', tool, '/rds/chr', chr, '_normal_sample.', tool, '.rds')
  
  if (to_parse == TRUE){
    if (tool == 'haplotypecaller'){
      rds_vcf = parse_HaplotypeCaller(file = vcf_file, out_file = rds_file, save = T)[[paste0(spn, '_normal_sample')]]$mutations
    } else if (tool == 'freebayes') {
      rds_vcf = parse_freebayes(file = vcf_file, out_file = rds_file, save = T, cutoff = 0.3)[[paste0(spn, '_normal_sample')]]$mutations
    }
  } else  {
    rds_vcf = readRDS(rds_file)[[paste0(spn, '_normal_sample')]]$mutations    
  }
  
  rds_vcf = rds_vcf %>% 
    dplyr::select(chr, from, to, ref, alt, NV, DP, BAF, QUAL, FILTER) %>% 
    dplyr::mutate(mutationID = paste(chr,from, sep = ':'))
  return(rds_vcf)
})

rds_vcf <- bind_rows(vcf_list)
rds_vcf_pass <- bind_rows(vcf_list) %>% dplyr::filter(FILTER == "PASS" & !is.na(FILTER))

# rds_ProCESS <- readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/', spn, '/races/seq_results_muts_merged_coverage_30x.rds'))  %>% 
#   dplyr::filter(classes == 'germinal') %>%
#   dplyr::mutate(chr = paste0('chr', chr)) %>% 
#   dplyr::filter(chr != 'chrX') %>% 
#   dplyr::filter(chr != 'chrY')  
# 
# ProCESS <- seq_to_long(rds_ProCESS) %>% 
#   dplyr::rename(BAF = VAF) %>% 
#   mutate(mutationID = paste(chr, from, sep = ':'))

ProCESS <- readRDS('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/ProCESS-examples/validation/Germline/ProCESS_germ.RDS')

sample_N <- 500000
p_baf_races <- ProCESS %>% 
  ungroup() %>% 
  sample_n(sample_N) %>% 
  ggplot(aes(x=BAF)) + 
  geom_histogram(binwidth=0.01) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  labs(x="BAF", 
       y="Count",
       title = "ProCESS BAF",
       subtitle = paste0(nrow(ProCESS), " mutations"),
       caption = paste0('sample ', sample_N, ' mutations'),
       color = "") 

p_dp_races <- ProCESS %>% 
  ungroup() %>% 
  sample_n(sample_N) %>% 
  ggplot(aes(x=DP)) + 
  geom_histogram(binwidth=1) +
  geom_vline(aes(xintercept = median(DP)), color = "indianred") +
  theme_bw() +
  labs(x="DP", 
       y="Count",
       title = "ProCESS coverage",
       subtitle = paste0(nrow(ProCESS), " mutations; Median coverage = ",median(ProCESS$DP)),
       caption = paste0('sample ', sample_N, ' mutations'),
       color = "") 

p_baf_caller <-  rds_vcf_pass %>% 
  sample_n(sample_N) %>% 
  ggplot(mapping = aes(x=BAF)) +
  geom_histogram(binwidth=0.01) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  labs(x="BAF", 
       y="Count",
       title = paste0(tool," BAF"),
       subtitle = paste0(nrow(rds_vcf_pass), " PASS mutations"),
       caption = paste0('sample ', sample_N, ' mutations'),
       color = "")

p_dp_caller <- rds_vcf_pass %>% 
  sample_n(sample_N) %>% 
  ggplot(mapping = aes(x=DP)) +
  geom_histogram(binwidth=1) +
  geom_vline(aes(xintercept = median(DP)), color = "indianred") +
  theme_bw() +
  labs(x="DP", 
       y="Count",
       title = paste0(tool ," coverage"),
       subtitle = paste0(nrow(rds_vcf_pass), " PASS mutations; Median coverage = ", median(rds_vcf$DP)),
       caption = paste0('sample ', sample_N, ' mutations'),
       color = "") 

merged_df <- merge_datasets(rds_vcf, ProCESS)
p_filter_dist = plot_filter_distribution(merged_df)

baf_differences = plot_baf_difference(merged_df)
cov_differences = plot_cov_difference(merged_df)

colors <- get_colors(merged_df)

#sample_N <- 1e5
merged_df_filter <- merged_df %>% ungroup() %>% sample_n(sample_N)

p_scatter_BAF_all = plot_scatter_with_corr(merged_df_filter, "BAF.races", "BAF.caller") +
  ggtitle("BAF correlation", subtitle = paste0(nrow(merged_df), " total mutations")) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "BAF races",
    y = paste0("BAF ", tool),
    caption = paste0('sample ', sample_N, ' mutations')
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p_scatter_BAF_all = ggExtra::ggMarginal(p_scatter_BAF_all, type = "boxplot", groupFill = TRUE, groupColour = TRUE)
p_scatter_BAF_all = ggplotify::as.ggplot(p_scatter_BAF_all)

p_scatter_DP_all = plot_scatter_with_corr(merged_df_filter, "DP.races", "DP.caller") +
  ggtitle("DP correlation", subtitle = paste0(nrow(merged_df), " total mutations")) +
  theme(legend.position = "bottom") +
  labs(
    x = "DP races",
    y = paste0("DP ", tool),
    caption = paste0('sample ', sample_N, ' mutations')
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p_scatter_DP_all = ggExtra::ggMarginal(
  p_scatter_DP_all,
  type = "boxplot",      # For density plots
  groupFill = TRUE, groupColour = TRUE
)
p_scatter_DP_all = ggplotify::as.ggplot(p_scatter_DP_all)

y_true = as.numeric(factor((merged_df$BAF.races > 0), levels=c(FALSE, TRUE))) - 1
y_pred = as.numeric(factor((merged_df$BAF.caller > 0), levels=c(FALSE, TRUE))) - 1

p_venn_all = plot_venn_diagram(merged_df, tool) +
  ggtitle("All called mutations")
metrics_all = compute_metrics(y_true, y_pred)

# PASS mutations
merged_df = merge_datasets(rds_vcf_pass, ProCESS)
merged_df_filter <- merged_df %>% ungroup() %>% sample_n(sample_N)

p_scatter_BAF_pass = plot_scatter_with_corr(merged_df_filter, "BAF.races", "BAF.caller") +
  ggtitle("BAF correlation", subtitle = paste0(nrow(merged_df), " total mutations")) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "BAF races",
    y = paste0("BAF ", tool),
    caption = paste0('sample ', sample_N, ' mutations')
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p_scatter_BAF_pass = ggExtra::ggMarginal(p_scatter_BAF_pass, type = "boxplot", groupFill = TRUE, groupColour = TRUE)
p_scatter_BAF_pass = ggplotify::as.ggplot(p_scatter_BAF_pass)

p_scatter_DP_pass = plot_scatter_with_corr(merged_df_filter, "DP.races", "DP.caller") +
  ggtitle("DP correlation", subtitle = paste0(nrow(merged_df), " total mutations")) +
  theme(legend.position = "bottom") +
  labs(
    x = "DP races",
    y = paste0("DP ", tool),
    caption = paste0('sample ', sample_N, ' mutations')
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p_scatter_DP_pass = ggExtra::ggMarginal(
  p_scatter_DP_pass,
  type = "boxplot",      # For density plots
  groupFill = TRUE, groupColour = TRUE
)
p_scatter_DP_pass = ggplotify::as.ggplot(p_scatter_DP_pass)

y_true = as.numeric(factor((merged_df$BAF.races > 0), levels=c(FALSE, TRUE))) - 1
y_pred = as.numeric(factor((merged_df$BAF.caller > 0), levels=c(FALSE, TRUE))) - 1

p_venn_pass = plot_venn_diagram(merged_df, tool) +
  ggtitle("Only PASS mutations")
metrics_pass = compute_metrics(y_true, y_pred)

metrics = dplyr::bind_rows(
  metrics_all %>% tidyr::pivot_longer(cols = colnames(metrics_all)) %>% dplyr::mutate(Mutations = "All"),
  metrics_pass %>% tidyr::pivot_longer(cols = colnames(metrics_pass)) %>% dplyr::mutate(Mutations = "Only Pass")  
)

p_metrics = metrics %>% 
  ggplot(mapping = aes(x=name, y=value, fill=Mutations)) +
  geom_col(position = "dodge") +
  theme_bw() +
  ylim(c(0,1))

design = "
AABBC
AABBC
DDEEF
DDEEF
GGHHI
GGHHI
LLMMN
LLMMN
OOOPP
"

title = paste0(spn, ", calls by ", tool)
report_plot = free(p_dp_races) + free(p_dp_caller) + free(cov_differences) +
  free(p_baf_races) + free(p_baf_caller) + free(baf_differences) +
  free(p_scatter_DP_all) + free(p_scatter_BAF_all) + free(p_venn_all) +
  free(p_scatter_DP_pass) + free(p_scatter_BAF_pass) + free(p_venn_pass) +
  free(p_metrics) + free(p_filter_dist) +
  plot_layout(design = design) +
  plot_annotation(title)

ggsave(plot = report_plot, filename = paste0(dir,'/', spn, '/normal_sample/', tool, '/plot/', tool, '_normal.png'), dpi = 400, width = 15, height = 20, units = 'in')

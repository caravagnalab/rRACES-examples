library(tidyverse)
library(vcfR)
library(optparse)
library(caret)

setwd('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/rRACES-examples/validation/Germline/')
source('vcf_parser.R')
source('utils.R')

option_list <- list( 
  make_option(c("-i", "--input"), type="character", default='/orfeo/LTS/LADE/LT_storage/lvaleriani/races/validation_data', help="path to input data"),
  make_option(c("-s", "--SPN"), type="character", default='SPN01', help="SPN name"),
  make_option(c("-t", "--tool"), type="character", default='freebayes', help="variant calling tool"),
  make_option(c("-C", "--chr"), type="character", default='22', help="chromosome number"),
  make_option(c("-o", "--output"), type="character", default='/orfeo/LTS/LADE/LT_storage/lvaleriani/races/validation_data', help="path to output directory")
)
to_parse = TRUE

param <- parse_args(OptionParser(option_list=option_list))
dir <- param$input
spn <- param$SPN
tool <- param$tool
output <- param$output
chr <- param$chr

vcf_file = paste0(dir,'/', spn, '/normal_sample/', tool, '/vcf/chr', chr, '_normal_sample.', tool, '.vcf.gz')
rds_file = paste0(dir,'/', spn, '/normal_sample/', tool, '/rds/chr', chr, '_normal_sample.', tool, '.rds')

if (to_parse == TRUE){
  if (tool == 'haplotypecaller'){
    rds_vcf <- parse_HaplotypeCaller(file = vcf_file, out_file = rds_file, save = T)[[paste0(spn, '_normal_sample')]]$mutations
  } else if (tool == 'freebayes') {
    rds_vcf <- parse_freebayes(file = vcf_file, out_file = rds_file, save = T, cutoff = 0.3)[[paste0(spn, '_normal_sample')]]$mutations
  }
} else  {
  rds_vcf <- readRDS(rds_file)[[paste0(spn, '_normal_sample')]]$mutations    
}

rds_vcf <- rds_vcf %>% select(chr, from, to, ref, alt, NV, DP, BAF, QUAL, FILTER) %>% mutate(mutationID = paste(chr,from, sep = ':'))

c <- chr
rds_rRACES <- readRDS(paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/', spn, '/races/seq_results_muts_merged_coverage_30x.rds')) %>% filter(chr == c)

rRACES <- seq_to_long(rds_rRACES) %>% 
  filter(classes == 'germinal') %>% 
  rename(BAF = VAF) %>% 
  mutate(chr = paste0('chr', chr)) %>% 
  mutate(mutationID = paste(chr,from, sep = ':'))

p_baf_races <- rRACES %>% 
  ggplot(aes(x=BAF)) + 
  geom_histogram(binwidth=0.01) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  labs(x="BAF", y="Count") +
  ggtitle("rRACES", subtitle = paste0(nrow(rRACES), " total mutations"))

p_baf_caller <-  rds_vcf %>% 
  dplyr::filter(FILTER == "PASS") %>% 
  ggplot(mapping = aes(x=BAF)) +
  geom_histogram(binwidth=0.01) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  labs(x="BAF", y="Count") +
  ggtitle(tool, subtitle = paste0(sum(rds_vcf$FILTER == "PASS"), " PASS mutations"))

# all
merged_df <- merge(rRACES, rds_vcf, by = "mutationID", all = TRUE, suffix=c(".races",".caller")) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .))

p_filter_dist = plot_filter_distribution(merged_df, log_scale = TRUE) + ggtitle("Distribution of caller flags")

p_scatter_BAF_all = plot_scatter_with_corr(merged_df, "BAF.races", "BAF.caller") +
  ggtitle("BAF correlation w.o filters", subtitle = paste0(nrow(merged_df), " total mutations")) +
  xlim(c(0,1)) +
  ylim(c(0,1))

p_scatter_DP_all = plot_scatter_with_corr(merged_df, "DP.races", "DP.caller") +
  ggtitle("DP correlation w.o filters", subtitle = paste0(nrow(merged_df), " total mutations"))

y_true = as.numeric(factor(!is.na(merged_df$BAF.races), levels=c(FALSE, TRUE))) - 1
y_pred = as.numeric(factor(!is.na(merged_df$BAF.caller), levels=c(FALSE, TRUE))) - 1
p_confusion_all = plot_confusion_matrix(y_true, y_pred)

metrics_all <- compute_metrics(y_true, y_pred) %>% 
  dplyr::mutate(Caller = tool, Filters = "Off")


# pass
merged_df = merge(rRACES, rds_vcf %>% dplyr::filter(FILTER=="PASS"), by = "mutationID", all = TRUE, suffix=c(".races",".caller")) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .))

p_scatter_BAF_pass = plot_scatter_with_corr(merged_df, "BAF.races", "BAF.caller") +
  ggtitle("BAF correlation w filters", subtitle = paste0(nrow(merged_df), " total mutations")) +
  xlim(c(0,1)) +
  ylim(c(0,1))

p_scatter_DP_pass = plot_scatter_with_corr(merged_df, "DP.races", "DP.caller") +
  ggtitle("DP correlation w filters", subtitle = paste0(nrow(merged_df), " total mutations"))

y_true = as.numeric(factor(!is.na(merged_df$BAF.races), levels=c(FALSE, TRUE))) - 1
y_pred = as.numeric(factor(!is.na(merged_df$BAF.caller), levels=c(FALSE, TRUE))) - 1
p_confusion_pass = plot_confusion_matrix(y_true, y_pred)

metrics_pass <- compute_metrics(y_true, y_pred) %>% 
  dplyr::mutate(Caller = tool, Filters = "On")

# Build report
metrics = dplyr::bind_rows(metrics_all, metrics_pass) %>% 
  tidyr::pivot_longer(!c(Filters, Caller))

p_metrics = metrics %>% 
  ggplot(mapping = aes(x=name, y=value, fill=Filters)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(x="", y="Value") +
  ggtitle("Metrics comparison")

design = "
AAABBB
CCCDDD
EEFFGG
HHIILL
"

report_plot = free(p_baf_races) + free(p_filter_dist) + 
  free(p_baf_caller) + free(p_metrics) + 
  free(p_scatter_BAF_all) + free(p_scatter_DP_all) + free(p_confusion_all) +
  free(p_scatter_BAF_pass) + free(p_scatter_DP_pass) + free(p_confusion_pass) +
  plot_layout(design = design) & 
  theme(legend.position = "bottom")  

ggsave(plot = report_plot, filename = paste0(dir,'/', spn, '/normal_sample/', tool, '/plot/', tool, '_normal_chr', chr, '.png'), dpi = 400, width = 10, height = 15, units = 'in')

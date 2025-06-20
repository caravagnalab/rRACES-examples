rm(list = ls())
options(bitmapType='cairo')
library(optparse)
library(tidyverse)

source('utils.R')

source('../../getters/process_getters.R')
source('../../getters/tumourevo_getters.R')

option_list <- list( 
  make_option(c("-v", "--variantcaller"), type="character", default='mutect2', help="variantcaller"),
  make_option(c("-c", "--cnacaller"), type="character", default='ascat', help="cna caller")
)

param <- parse_args(OptionParser(option_list=option_list))

SPN <- paste0('SPN0', seq(1,4))
COV <- c(50, 100, 150, 200)
PUR <- c(0.3, 0.6, 0.9)

full_table <- tibble()
for (spn in SPN){
  out <- paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/', spn, '/validation/tumourevo/CNAqc/')
  dir.create(out, showWarnings = F, recursive = T)
  print(spn)
  
  samples <- get_sample_names(spn)
  for (cov in COV){
    print(cov)
    for (pur in PUR){
      print(pur)
    
      cnaqc_obj <- lapply(samples, FUN = function(sample){
        file <- get_tumourevo_qc(spn = spn, coverage = cov, purity = pur, tool = 'CNAqc', vcf_caller = param$variantcaller, cna_caller = param$cnacaller, sample = sample)
        if (length(file) > 0){
          if (file.exists(file$qc_rds)){
          readRDS(file$qc_rds)
          }
        }
      })
      names(cnaqc_obj) <- samples
      
      if (!is.null(cnaqc_obj[[1]])){
        stats_table <- get_statistics_qc(cnaqc_obj, purity = pur)
        table <- stats_table$table %>% tibble() %>% select(-gg) %>% filter(info != 'sample') %>% tidyr::pivot_wider(values_from = values, names_from = info, names_repair = 'minimal')
        full_table <- bind_rows(full_table, table %>% mutate(spn = spn, coverage = cov, vcf_caller = param$variantcaller, cna_caller = param$cnacaller))
      }
      #plot_simple <- plot_qc(cnaqc_obj, type = 'simple_clonal')
      #plot_complex <- plot_qc(cnaqc_obj, type = 'complex_clonal')
    }
  }
}

qc <- ggplot() + 
  ggplot2::geom_tile(
    data = full_table,
    ggplot2::aes(
      y = '',
      x = sample,
      fill = QC,
      width = .8,
      height = .8
    )
  ) +
  ggh4x::facet_nested(coverage+true_purity ~ spn, scales = 'free_x')  +
  ggplot2::scale_fill_manual(values = c(`FAIL` = 'indianred3', `PASS` = 'seagreen', `NA` = 'gainsboro'))  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('')



df_clean <- full_table %>%
  mutate(
    sample_id = sample,
    n_clonal_cnas = as.numeric(n_clonal_cnas),
    n_subclonal_cnas = as.numeric(n_subclonal_cnas),
    passed_clonal_cnas = as.numeric(passed_clonal_cnas0),
  ) %>%
  mutate(
    total_cnas = n_clonal_cnas,
    passed_cnas = passed_clonal_cnas,
    failed_cnas = total_cnas - passed_cnas,
    pct_passed = passed_cnas / total_cnas * 100,
    pct_failed = 100 - pct_passed
  )

df_long <- df_clean %>%
  select(sample_id, spn, pct_passed, pct_failed, true_purity, coverage) %>%
  pivot_longer(cols = starts_with("pct_"), names_to = "status", values_to = "percentage") %>%
  mutate(status = recode(status, pct_passed = "Passed", pct_failed = "Failed"))


cnas <- ggplot(df_long, aes(x = sample_id, y = percentage, fill = status)) +
  geom_bar(stat = "identity") +
  ggh4x::facet_nested(coverage+true_purity ~ spn, scales = 'free_x')  +
  labs(x = "sample", y = "Percentage of segment", fill = "Segment status") +
  scale_fill_manual(values = c("Passed" = "seagreen", "Failed" = "gainsboro")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plt <- qc + cnas + plot_layout(nrow = 2)
ggsave(filename = 'validation_CNAqc.png', plot = plt, width = 7, height = 7, units = 'in', dpi = 600)

rm(list = ls())
options(bitmapType='cairo')
library(tidyverse)
library(optparse)
library(caret)
library(dplyr)
library(patchwork)

option_list <- list( 
  make_option(c("-v", "--variantcaller"), type="character", default='mutect2', help="variantcaller")
)
param <- parse_args(OptionParser(option_list=option_list))

SPN <- paste0('SPN0', seq(1,7))
COV <- c(50, 100, 150, 200)
PUR <- c(0.3, 0.6, 0.9)

final_table <- tibble()
for (spn in SPN){
  for (cov in COV){
    for (pur in PUR){
      base <- paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/',spn, '/tumourevo/',cov,'x_',pur,'p_',param$variantcaller,'_ascat/QC/tinc/SCOUT/',spn, '/')
      if (file.exists(base)){
        print(base)
        
        samples <- list.dirs(base, full.names = F, recursive = F)
        table <- lapply(samples, FUN = function(sample){
          files <- list.files(paste0(base, sample))
          data <- readRDS(paste0(base, sample, '/', files[grepl('fit', files)]))
          s <- str_split(sample, '_') %>% unlist()
          tmp <- tibble(TIT = data$TIT, TIN = data$TIN, sample = paste0(s[2], '_', s[3]))
          return(tmp)
        }) %>% bind_rows()
        
        table$coverage = as.numeric(cov)
        table$purity = as.numeric(pur)
        table$vc = param$variantcaller
        table$cnc = 'ascat'
        table$SPN = spn
        table$N = length(samples)
        table <- table %>% mutate(error = abs(purity-TIT))
        #saveRDS(object = table, file = paste0(out, param$cov,'x_',param$pur,'p_',param$variantcaller,'_ascat.rds'))
        
        final_table <- bind_rows(final_table, table)
      }
    }
  }
}

tit <- final_table %>% 
  ggplot() +
  geom_abline(linewidth = 0.5, col = 'gray') +
  geom_jitter(aes(x = TIT, y = purity, col = SPN), size = 3) +
  geom_point(aes(x = purity, y = purity), size = 3, shape = 8) +
  facet_grid(coverage ~ purity) + 
  xlab('TIT (TINC)') +
  ylab('purity (ProCESS)') + 
  scale_color_manual(values = c('steelblue', 'seagreen', 'goldenrod')) +
  xlim(0,1) +
  ylim(0,1) +
  theme_minimal()

ggsave(filename = 'validation_TINC.png', plot = tit, width = 7, height = 3, units = 'in', dpi = 600)

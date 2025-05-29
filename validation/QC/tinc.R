rm(list = ls())
options(bitmapType='cairo')
library(tidyverse)
library(optparse)
library(caret)
library(dplyr)
library(patchwork)

option_list <- list( 
  make_option(c("-s", "--SPN"), type="character", default='SPN03', help="SPN name"),
  make_option(c("-c", "--cov"), type="character", default='50', help="coverage"),
  make_option(c("-p", "--pur"), type="character", default='0.6', help="purity"),
  make_option(c("-v", "--variantcaller"), type="character", default='mutect2', help="variantcaller")
)

param <- parse_args(OptionParser(option_list=option_list))
base <- paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/',param$SPN, '/tumourevo/',param$cov,'x_',param$pur,'p_',param$variantcaller,'_ascat/QC/tinc/SCOUT/', param$SPN, '/')
samples <- list.dirs(base, full.names = F, recursive = F)

out <- paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/', param$SPN, '/validation/tumourevo/tinc/')
dir.create(out, showWarnings = F, recursive = T)

table <- lapply(samples, FUN = function(sample){
  files <- list.files(paste0(base, sample))
  data <- readRDS(paste0(base, sample, '/', files[grepl('fit', files)]))
  s <- str_split(sample, '_') %>% unlist()
  tmp <- tibble(TIT = data$TIT, TIN = data$TIN, sample = paste0(s[2], '_', s[3]))
  table <- bind_rows(table, tmp)
}) %>% bind_rows()

table$coverage = as.numeric(param$cov)
table$purity = as.numeric(param$pur)
table$vc = param$variantcaller
table$cnc = 'ascat'
table <- table %>% mutate(error = abs(purity-TIT))
saveRDS(object = table, file = paste0(out, param$cov,'x_',param$pur,'p_',param$variantcaller,'_ascat.rds'))

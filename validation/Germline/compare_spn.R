rm(list = ls())
options(bitmapType='cairo')
library(tidyverse)
library(vcfR)
library(optparse)
library(caret)
library(dplyr)
library(patchwork)

source('utils.R')
source('../../getters/sarek_getters.R')
source('../../getters/process_getters.R')

dir <- '/orfeo/scratch/cdslab/shared/SCOUT/'

all_metric <- tibble()
all_baf <- tibble()

spn <- paste0('SPN0', 1:7)
for (s in spn){

  input <- paste0(dir, s, "/validation/germline/report")
  if (file.exists(file.path(input, 'normal_metrics.rds'))){
    
    rds <- readRDS(file.path(input, 'normal_metrics.rds'))
    metric <- rds$report_metrics
    baf <- rds$baf_metric

    all_metric <- bind_rows(all_metric, metric) 
    all_baf <- bind_rows(all_baf, baf) 
    
  }
}
# do some plots
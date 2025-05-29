rm(list = ls())
options(bitmapType='cairo')
library(tidyverse)
library(optparse)
library(CNAqc)
library(dplyr)
library(patchwork)

source('utils.R')

option_list <- list( 
  make_option(c("-s", "--SPN"), type="character", default='SPN01', help="SPN name"),
  make_option(c("-c", "--cov"), type="character", default='50', help="coverage"),
  make_option(c("-p", "--pur"), type="character", default='0.6', help="purity"),
  make_option(c("-v", "--variantcaller"), type="character", default='mutect2', help="variantcaller")
)

param <- parse_args(OptionParser(option_list=option_list))
base <- paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/',param$SPN, '/tumourevo/',param$cov,'x_',param$pur,'p_',param$variantcaller,'_ascat/QC/CNAqc/SCOUT/', param$SPN, '/')
samples <- list.dirs(base, full.names = F, recursive = F)

out <- paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/', param$SPN, '/validation/tumourevo/CNAqc/')
dir.create(out, showWarnings = F, recursive = T)

cnaqc_obj <- lapply(samples, FUN = function(sample){
  files <- list.files(paste0(base, sample))
  readRDS(paste0(base, sample, '/', files[grepl('qc.rds', files)]))
})
names(cnaqc_obj) <- samples

# create res
stats_table <- get_statistics_qc(cnaqc_obj, purity = param$pur)
plot_simple <- plot_qc(cnaqc_obj, type = 'simple_clonal')
plot_complex <- plot_qc(cnaqc_obj, type = 'complex_clonal')

rm(list = ls())

setwd("/orfeo/cephfs/scratch/cdslab/gsantacatterina/ProCESS-examples/validation/Somatic")

require(tidyverse)
library(optparse)
source("utils/mutect_utils.R")
source("utils/process_utils.R")
source("utils/freeBayes_utils.R")
source("utils/strelka_utils.R")
source("../../getters/process_getters.R")
source("../../getters/sarek_getters.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN01'),
                    make_option(c("--purity"), type = "character", default = '0.6'),
                    make_option(c("--coverage"), type = "character", default = '100'),
                    make_option(c("--chr"), type = "character", default = '22')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT'
spn_id = opt$spn_id
coverage = opt$coverage
purity = opt$purity
chromosome <- opt$chr
outdir <-  file.path(data_dir,spn_id,"/validation/somatic/")

gt_path = get_mutations(spn = spn_id, 
                        base_path = data_dir, 
                        coverage = coverage, 
                        purity = purity, 
                        type = "tumour")

message("PARSING RACES")
process_seq_results(gt_path, spn_id, purity, coverage, chromosome, base_path = data_dir, outdir = outdir)

message("PARSING MUTECT2")
process_mutect2_results(gt_path, spn_id, purity, coverage, chromosome, base_path = data_dir, outdir = outdir)

message("PARSING FREEBAYES")
process_freebayes_results(gt_path, spn_id, purity, coverage, chromosome, base_path = data_dir, outdir = outdir, pass_quality = 20, min_vaf = 0.01, max_normal_vaf = 0.02)

message("PARSING STRELKA")
process_strelka_results(gt_path, spn_id, purity, coverage, chromosome, base_path = data_dir, outdir = outdir)

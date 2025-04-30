rm(list = ls())
require(tidyverse)
source("utils/mutect_utils.R")
source("utils/races_utils.R")
source("utils/freeBayes_utils.R")
source("utils/strelka_utils.R")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript process_mutect2_res.R <chromosome>")
}

# Get the chromosome from the command line arguments
chromosome <- args[1]

gt_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/races/purity_0.3/seq_results_muts_merged_coverage_50x.rds"

message("PARSING RACES")
process_seq_results(gt_path, chromosome)
message("PARSING MUTECT2")
process_mutect2_results(gt_path, chromosome)
message("PARSING FREEBAYES")
process_freebayes_results(gt_path, chromosome, pass_quality = 20, min_vaf = 0.01, max_normal_vaf = 0.02)
message("PARSING STRELKA")
process_strelka_results(gt_path, chromosome)

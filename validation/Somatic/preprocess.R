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

path_to_seq = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/sequencing/tumour/purity_0.3/data/mutations/seq_results_muts_merged_coverage_50x.rds"
mutect_vcfs_dir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/sarek/50x_0.3p/variant_calling/mutect2/SPN01/"
freebayes_vcfs_dir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/sarek/50x_0.3p/variant_calling/freebayes/"
strelka_vcfs_dir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/sarek/50x_0.3p/variant_calling/strelka/"
outdir = "/orfeo/cephfs/scratch/cdslab/gsantacatterina/rRACES_test/outdir"

message("PARSING RACES")
process_seq_results(path_to_seq, chromosome, outdir)
message("PARSING MUTECT2")
process_mutect2_results(path_to_seq, chromosome, outdir, mutect_vcfs_dir)
message("PARSING FREEBAYES")
process_freebayes_results(path_to_seq, chromosome, outdir, freebayes_vcfs_dir, pass_quality = 20, min_vaf = 0.01, max_normal_vaf = 0.02)
message("PARSING STRELKA")
process_strelka_results(path_to_seq, chromosome, outdir, strelka_vcfs_dir)

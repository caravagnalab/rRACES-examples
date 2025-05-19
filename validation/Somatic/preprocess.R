rm(list = ls())
require(tidyverse)
library(optparse)
source("utils/mutect_utils.R")
source("utils/process_utils.R")
source("utils/freeBayes_utils.R")
source("utils/strelka_utils.R")



option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN01'),
                    make_option(c("--purity"), type = "character", default = '0.6'),
                    make_option(c("--coverage"), type = "character", default = '100'),
                    make_option(c("--chr"), type = "character", default = '22')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'
spn_id = opt$spn_id
coverage = opt$coverage
purity = opt$purity
chromosome <- opt$chr


path_to_seq <- paste0(data_dir,spn_id,"/sequencing/tumour/purity_",purity,"/data/mutations/seq_results_muts_merged_coverage_",coverage,"x.rds")
mutect_vcfs_dir <- paste0(data_dir,spn_id,"/sarek/",coverage,"x_",purity,"p","/variant_calling/mutect2/",spn_id)
freebayes_vcfs_dir <- paste0(data_dir,spn_id,"/sarek/",coverage,"x_",purity,"p","/variant_calling/freebayes/")
strelka_vcfs_dir <- paste0(data_dir,spn_id,"/sarek/",coverage,"x_",purity,"p","/variant_calling/strelka/")
outdir <-  paste0(data_dir,spn_id,"/validation/somatic/")

message("PARSING RACES")
process_seq_results(path_to_seq, chromosome, outdir)
chromosome <- paste0("chr",chromosome)
message("PARSING MUTECT2")
process_mutect2_results(path_to_seq, chromosome, outdir, mutect_vcfs_dir)
message("PARSING FREEBAYES")
process_freebayes_results(path_to_seq, chromosome, outdir, freebayes_vcfs_dir, pass_quality = 20, min_vaf = 0.01, max_normal_vaf = 0.02)
message("PARSING STRELKA")
process_strelka_results(path_to_seq, chromosome, outdir, strelka_vcfs_dir)

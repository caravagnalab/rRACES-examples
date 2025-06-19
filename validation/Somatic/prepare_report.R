
rm(list = ls())
options(bitmapType='cairo')
require(tidyverse)
library(optparse)

source("../../getters/process_getters.R")
source("utils/utils.R")
source("utils/plot_utils.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN04'),
                    make_option(c("--purity"), type = "character", default = '0.3'),
                    make_option(c("--coverage"), type = "character", default = '50'))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'
spn_id = opt$spn_id
coverage = opt$coverage
purity = opt$purity

# gender_file <- read.table(file = paste0(data_dir,spn_id,"/process/subject_gender.txt"),header = FALSE,col.names = "gender")
# gender <- gender_file$gender
gender = get_process_gender(spn = spn_id)
if (gender=="XX"){
  chromosomes = c(paste0('chr',1:22), 'chrX')
} else {
  chromosomes = c(paste0('chr',1:22), 'chrX', 'chrY')
}

# INPUT PARAMATERS ####
callers = c("mutect2", "strelka", "freebayes")
# chromosomes = paste0("chr", c(1:22, "X", "Y"))
min_vaf = .02
mut_types = c("INDEL", "SNV")
comb = list(PI = purity, COV = coverage)
#path_to_cna <- paste0(data_dir,spn_id,"/process/cna_data/")
#samples = gsub(pattern="_cna.rds",replacement="",x=list.files(path_to_cna,pattern="_cna.rds"))
samples = get_sample_names(spn = spn_id)

# path_to_seq = get_mutations(spn = spn_id, 
#                         base_path = data_dir, 
#                         coverage = coverage, 
#                         purity = purity, 
#                         type = "tumour")

#path_to_seq <- paste0(data_dir,spn_id,"/sequencing/tumour/purity_",purity,"/data/mutations/seq_results_muts_merged_coverage_",coverage,"x.rds")

input_dir <-  paste0(data_dir,spn_id,"/validation/somatic/")
outdir <- paste0(data_dir,spn_id,"/validation/somatic/report")
dir.create(outdir, recursive = T, showWarnings = F)

# Preparing report
message("Parsing combination: purity=", purity, ", cov=", coverage)
caller = "mutect2"
sample_id = samples[1]
mut_type = "SNV"
for (caller in callers) {
  message("  Using caller: ", caller)
  for (sample_id in samples) {
    message("    Working with sample : ", sample_id)
    for (mut_type in mut_types) {
      message("      Considering only ", mut_type)
      
      # spn <- gsub(".*SCOUT/(SPN[0-9]+).*", "\\1", path_to_seq)
      # purity <- gsub(".*purity_([0-9.]+).*", "\\1", path_to_seq)
      # coverage <- gsub(".*coverage_([0-9]+).*", "\\1", path_to_seq)
      combination = paste0(coverage, "x_", purity, "p")
      process_folder_path <- file.path(input_dir, spn_id, combination, "process", sample_id, mut_type)
      caller_folder_path = file.path(input_dir, spn_id, combination, caller, sample_id, mut_type)
      
      # Get ground truth
      gt_res = lapply(chromosomes, function(chromosome){
        gt_path = file.path(process_folder_path, paste0(chromosome,".rds"))
        readRDS(gt_path)
      }) %>% do.call("bind_rows", .)
      
      # Get caller res
      caller_res = lapply(chromosomes, function(chromosome) {
        caller_path = file.path(caller_folder_path, paste0(chromosome,".rds"))
        readRDS(caller_path)  
      }) %>% do.call("bind_rows", .) %>% 
        dplyr::filter(!is.na(VAF))
      
      sample_info = list(caller_name=caller, sample_id=sample_id, mut_type=mut_type, spn=spn_id, purity=purity, coverage=coverage)
      report = get_report(seq_res_long = gt_res, 
                          caller_res = caller_res, 
                          sample_info = sample_info, 
                          min_vaf = min_vaf)
      
      report_path = file.path(caller_folder_path, "report.png")
      metrics_path = file.path(caller_folder_path, "metrics.rds")
      ggsave(report_path, plot = report$report_plot, width = 15, height = 20, units = "in", dpi = 400)
      saveRDS(list(report_metrics=report$report_metrics, vaf_comparison=report$vaf_comparison), metrics_path)
      
      filename = paste(spn_id, combination, caller, sample_id, mut_type, sep = '_')
      file_path = file.path(outdir, filename)
      ggsave(paste0(file_path, '.png'), plot = report$report_plot, width = 15, height = 20, units = "in", dpi = 400)
      
      message("        Report done!")
    }
  }
}

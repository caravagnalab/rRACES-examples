
rm(list = ls())
require(tidyverse)
source("utils/utils.R")
source("utils/plot_utils.R")

# INPUT PARAMATERS ####
callers = c("mutect2", "strelka", "freebayes")
chromosomes = paste0("chr", c(1:22, "X", "Y"))
min_vaf = .02
purity = 0.3
coverage = 100
mut_types = c("INDEL", "SNV")
comb = list(PI = purity, COV = coverage)
spn = "SPN01"
samples = c("SPN01_1.1", "SPN01_1.2", "SPN01_1.3")

path_to_seq = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/sequencing/tumour/purity_0.3/data/mutations/seq_results_muts_merged_coverage_50x.rds"
input_dir = "/orfeo/cephfs/scratch/cdslab/gsantacatterina/rRACES_test/outdir"
outdir = "/orfeo/cephfs/scratch/cdslab/gsantacatterina/rRACES_test/report_dir"
dir.create(outdir, recursive = T)

# Preparing report
message("Parsing combination: purity=", purity, ", cov=", coverage)

for (caller in callers) {
  message("  Using caller: ", caller)
  for (sample_id in samples) {
    message("    Working with sample : ", sample_id)
    
    for (mut_type in mut_types) {
      message("      Considering only ", mut_type)
      
      spn <- gsub(".*SCOUT/(SPN[0-9]+).*", "\\1", path_to_seq)
      purity <- gsub(".*purity_([0-9.]+).*", "\\1", path_to_seq)
      coverage <- gsub(".*coverage_([0-9]+).*", "\\1", path_to_seq)
      combination = paste0(coverage, "x_", purity, "p")
      process_folder_path <- file.path(input_dir, spn, combination, "process", sample_id, mut_type)
      caller_folder_path = file.path(input_dir, spn, combination, caller, sample_id, mut_type)
      
      # Get ground truth
      gt_res = lapply(chromosomes, function(chromosome){
        #print(chromosome)
        gt_path = file.path(process_folder_path, paste0(chromosome,".rds"))
        readRDS(gt_path)
      }) %>% do.call("bind_rows", .)
      
      # Get caller res
      caller_res = lapply(chromosomes, function(chromosome) {
        #print(chromosome)
        caller_path = file.path(caller_folder_path, paste0(chromosome,".rds"))
        readRDS(caller_path)  
      }) %>% do.call("bind_rows", .) %>% 
        dplyr::filter(!is.na(VAF))
      
      sample_info = list(caller_name=caller, sample_id=sample_id, mut_type=mut_type, spn=spn, purity=purity, coverage=coverage)
      report = get_report(gt_res, caller_res, sample_info, min_vaf)
      
      report_path = file.path(caller_folder_path, "report.png")
      metrics_path = file.path(caller_folder_path, "metrics.rds")
      #report_path = paste0(spn,"/purity_",purity,"_coverage_",coverage,"x/",caller,"/",sample_id,"/",mut_type,"/report.png")
      #metrics_path = paste0(spn,"/purity_",purity,"_coverage_",coverage,"x/",caller,"/",sample_id,"/",mut_type,"/metrics.rds")
      
      ggsave(report_path, plot = report$report_plot, width = 15, height = 20, units = "in", dpi = 400)
      saveRDS(list(report_metrics=report$report_metrics, vaf_comparison=report$vaf_comparison), metrics_path)
      
      message("        Report done!")
    }
  }
}


rm(list = ls())
require(tidyverse)
source("utils/utils.R")
source("utils/plot_utils.R")

# sample_id = "SPN01_1.2"
# caller = "strelka"
# mut_type = "INDEL"
# purity = 0.3
# coverage = 50
chromosomes = paste0("chr", c(1:22, "X", "Y"))

callers = c("mutect2", "strelka", "freebayes")
min_vaf = .02
mut_types = c("INDEL", "SNV")

combinations = dplyr::tibble(
  COV = c(100, 50),
  PI = c(.6, .3)
)

spn = "SPN01"
samples = c("SPN01_1.1", "SPN01_1.2", "SPN01_1.3")

for (idx in 1:nrow(combinations)) {
  comb = combinations[idx,]
  purity = comb$PI
  coverage = comb$COV
  message("Parsing combination: purity=", purity, ", cov=", coverage)
  
  for (caller in callers) {
    message("  Using caller: ", caller)
    for (sample_id in samples) {
      message("    Working with sample : ", sample_id)
      
      for (mut_type in mut_types) {
        message("      Considering only ", mut_type)
        
        # Get ground truth
        gt_res = lapply(chromosomes, function(chromosome){
          #print(chromosome)
          gt_path = paste0(spn,"/purity_",purity,"_coverage_",coverage,"x/races/",sample_id,"/",mut_type,"/",chromosome,".rds")  
          readRDS(gt_path)
        }) %>% do.call("bind_rows", .)
        
        # Get caller res
        caller_res = lapply(chromosomes, function(chromosome) {
          #print(chromosome)
          caller_path = paste0(spn,"/purity_",purity,"_coverage_",coverage,"x/",caller,"/",sample_id,"/",mut_type,"/",chromosome,".rds")
          readRDS(caller_path)  
        }) %>% do.call("bind_rows", .) %>% 
          dplyr::filter(!is.na(VAF))
        
        sample_info = list(caller_name=caller, sample_id=sample_id, mut_type=mut_type, spn=spn, purity=purity, coverage=coverage)
        report = get_report(gt_res, caller_res, sample_info, min_vaf)
        
        report_path = paste0(spn,"/purity_",purity,"_coverage_",coverage,"x/",caller,"/",sample_id,"/",mut_type,"/report.png")
        metrics_path = paste0(spn,"/purity_",purity,"_coverage_",coverage,"x/",caller,"/",sample_id,"/",mut_type,"/metrics.rds")
        
        ggsave(report_path, plot = report$report_plot, width = 15, height = 20, units = "in", dpi = 400)
        saveRDS(list(report_metrics=report$report_metrics, vaf_comparison=report$vaf_comparison), metrics_path)
        
        message("        Report done!")
      }
    }
  }
}

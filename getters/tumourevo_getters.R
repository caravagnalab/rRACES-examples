
source("utils.R")

################################################################################
# DRIVERS
################################################################################
get_tumourevo_driver <- function(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    sample, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  MAIN_PATH <- dir_getter(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    base_path
  )
  
  MAIN_PATH <- file.path(MAIN_PATH, "driver_annotation/annotate_driver/SCOUT", spn)
  
  output <- list.files(MAIN_PATH, pattern = paste0("\\", sample, "_driver.rds$"), full.names = T, recursive=T)
  
  return(output)
}


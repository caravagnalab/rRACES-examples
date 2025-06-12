
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



################################################################################
# SUBCLONAL
################################################################################
get_tumourevo_subclonal <- function(
    spn, 
    coverage, 
    purity, 
    tool, 
    vcf_caller, 
    cna_caller, 
    sample, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  # quality control
  tool_list <- c("mobster", "pyclone", "ctree", "viber")
  if (!(tool %in% tool_list)) {
    stop("ERROR: wrong tool name!")
  }
  
  MAIN_PATH <- dir_getter(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    base_path
  )
  
  MAIN_PATH <- file.path(MAIN_PATH, "subclonal_deconvolution", tool, "SCOUT", spn)
  
  if (tool=="ctree") {
    MAIN_PATH <- file.path(MAIN_PATH, paste0(spn, "_", sample))
  } else if (tool=="mobster") {
    MAIN_PATH <- file.path(MAIN_PATH, paste0(spn, "_", sample))
  }
  
  output <- get_named_file_list(MAIN_PATH)
  
  return(output)
}



################################################################################
# QC
################################################################################
get_tumourevo_qc <- function(
    spn, 
    coverage, 
    purity, 
    tool, 
    vcf_caller, 
    cna_caller, 
    sample, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  # quality control
  tool_list <- c("CNAqc", "join_CNAqc", "tinc")
  if (!(tool %in% tool_list)) {
    stop("ERROR: wrong tool name!")
  }
  
  MAIN_PATH <- dir_getter(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    base_path
  )
  
  MAIN_PATH <- file.path(MAIN_PATH, "QC", tool, "SCOUT", spn)
  
  all_entries <- list.dirs(MAIN_PATH, full.names = FALSE, recursive = FALSE)
  matching_dir <- all_entries[grepl(sample, all_entries)]
  MAIN_PATH <- file.path(MAIN_PATH, matching_dir)
  output <- get_named_file_list(MAIN_PATH)
  
  return(output)
  
}



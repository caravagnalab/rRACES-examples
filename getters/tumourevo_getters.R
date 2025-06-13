
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
    output <- get_named_file_list(MAIN_PATH)
    output <- ctree_named_list(output)
  } else if (tool=="mobster") {
    MAIN_PATH <- file.path(MAIN_PATH, paste0(spn, "_", sample))
    output <- get_named_file_list(MAIN_PATH)
    output <- mobster_named_list(output)
  } else if (tool=="pyclone") {
    output <- get_named_file_list(MAIN_PATH)
    output <- pyclone_input_list(output)
  } else if (tool=="viber") {
    output <- get_named_file_list(MAIN_PATH)
    output <- viber_input_list(output)
  }
  return(output)
}


#TODO:
# - join_CNAqc : no output
# - CNAqc : format the output
# - tinc : format the output
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




################################################################################
# SIGNATURES
################################################################################
get_tumourevo_signatures <- function(
    spn, 
    coverage, 
    purity, 
    tool, 
    context, 
    vcf_caller, 
    cna_caller, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  # quality control
  tool_list <- c("SigProfiler", "SparseSignatures")
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
  
  MAIN_PATH <- file.path(MAIN_PATH, "signature_deconvolution")
  
  
  
  if (tool == "SigProfiler") {
    
    MAIN_PATH <- file.path(MAIN_PATH, tool, "SCOUT", "results", context, "Suggested_Solution")
    
    COSMIC <- file.path(MAIN_PATH, paste0("COSMIC_", context, "_Decomposed_Solution"))
    COSMIC_exposure <- paste0("Activities/COSMIC_", context, "_Activities.txt")
    COSMIC_signatures <- paste0("Signatures/COSMIC_", context, "_Signatures.txt")
    
    denovo <- file.path(MAIN_PATH, paste0(context, "_De-Novo_Solution"))
    denovo_exposure <- paste0("Activities/", context, "_De-Novo_Activities_refit.txt")
    denovo_signatures <- paste0("Signatures/", context, "_De-Novo_Signatures.txt")
    
    return(list(
      COSMIC_exposure=COSMIC_exposure, 
      COSMIC_signatures=COSMIC_signatures, 
      denovo_exposure=denovo_exposure, 
      denovo_signatures=denovo_signatures
    ))
    
  } else if (tool == "SparseSignatures") {
    
    MAIN_PATH <- file.path(MAIN_PATH, tool, "SCOUT")
    output <- file.path(MAIN_PATH, "SCOUT_nmf_Lasso_out.rds")
    return(output)
  } else {
    stop("ERROR: wrong tool name!")
  }
  
}


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
    
    mut_counts=paste0(MAIN_PATH, "/", tool, "/SCOUT/results/", context, "/Samples.txt")
    MAIN_PATH <- file.path(MAIN_PATH, tool, "SCOUT", "results", context, "Suggested_Solution")
    
    COSMIC <- file.path(MAIN_PATH, paste0("COSMIC_", context, "_Decomposed_Solution"))
    COSMIC_exposure <- paste0(COSMIC, "/Activities/COSMIC_", context, "_Activities.txt")
    COSMIC_signatures <- paste0(COSMIC, "/Signatures/COSMIC_", context, "_Signatures.txt")
    
    denovo <- file.path(MAIN_PATH, paste0(context, "_De-Novo_Solution"))
    denovo_exposure <- paste0(denovo, "/Activities/", context, "_De-Novo_Activities_refit.txt")
    denovo_signatures <- paste0(denovo, "/Signatures/", context, "_De-Novo_Signatures.txt")
    
    return(list(
      mut_counts=mut_counts,
      COSMIC_exposure=COSMIC_exposure, 
      COSMIC_signatures=COSMIC_signatures, 
      denovo_exposure=denovo_exposure, 
      denovo_signatures=denovo_signatures
    ))
    
  } else if (tool == "SparseSignatures") {
    
    MAIN_PATH <- file.path(MAIN_PATH, tool, "SCOUT")
    output <- list() 
    output[['nmf_Lasso_out']] <- file.path(MAIN_PATH, "SCOUT_nmf_Lasso_out.rds")
    output[['cv_means_mse']] <- file.path(MAIN_PATH, "SCOUT_cv_means_mse.rds")
    output[['best_params_config']] <- file.path(MAIN_PATH, "SCOUT_best_params_config.rds")
    
    return(output)
  } else {
    stop("ERROR: wrong tool name!")
  }
  
}


samples_table <- function(snapshot, sample_forest) {
  sim <- ProCESS::recover_simulation(snapshot)
  info = sim$get_samples_info() ## requested from either the simulation recovery or as saved table
  
  nodes = sample_forest$get_nodes()
  clones = nodes %>% 
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    dplyr::pull(mutant) %>% 
    unique()
  clones_of_origin = nodes %>%
    dplyr::filter(!is.na(sample)) %>% 
    dplyr::group_by(sample, mutant) %>% 
    dplyr::count(mutant) %>% 
    tidyr::pivot_wider(values_from = n, names_from = mutant, values_fill = 0) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(Sample_Type = sum(c_across(clones) == 0)) %>% 
    dplyr::mutate(Sample_Type = ifelse(Sample_Type == (length(clones)-1), "Monoclonal", "Polyclonal")) %>% 
    dplyr::mutate(Total_Cells = sum(c_across(clones), na.rm = T)) %>% 
    dplyr::mutate(across(all_of(clones), ~ round(.x/Total_Cells,2), .names = "{.col} proportion"))
  
  samples_tb = dplyr::full_join(info, clones_of_origin, by = join_by("name" == "sample")) %>% 
    dplyr::select(!c("xmin","xmax","ymin","ymax","id","tumour_cells","tumour_cells_in_bbox")) %>% 
    dplyr::rename(Sample_ID=name) %>% 
    dplyr::rename(Samping_Time=time) %>% 
    dplyr::mutate(Samping_Time=round(Samping_Time,2)) %>% 
    dplyr::arrange(Sample_ID)
  
  return(samples_tb)
}


get_sbs_exposures <- function(spn, base_path) {

  phylo_forest <- get_phylo_forest(spn = spn, base_path = base_path)
  sample_forest <- get_sample_forest(spn = spn, base_path = base_path)
  snapshot_path <- file.path(base_path, spn, "process", spn)

  # Load phyloforest object
  phylo_forest <- ProCESS::load_phylogenetic_forest(phylo_forest)
  # Load sample forest object
  samples_forest <- ProCESS::load_samples_forest(sample_forest)

  # Generate sample-level data
  samples_data <- samples_table(snapshot = snapshot_path, sample_forest = samples_forest)

  sample_ids <- samples_data %>%
    dplyr::pull(Sample_ID)

  # Get exposures
  exposure_prop <- phylo_forest$get_exposures()

  # Filter out time = 0 exposures if others exist
  if (any(exposure_prop$time != 0)) {
    exposure_prop_filtered <- exposure_prop %>% dplyr::filter(time != 0)
  } else {
    exposure_prop_filtered <- exposure_prop
  }

  # Match sample IDs to timepoints
  assign_sample_id <- function(t) {
    if (t == 0) {
      return(sample_ids)
    } else {
      match_idx <- which.min(abs(samples_data$Samping_Time - t))
      return(samples_data$Sample_ID[match_idx])
    }
  }

  # Expand each exposure row to a sample ID
  expanded_rows <- do.call(rbind, lapply(1:nrow(exposure_prop_filtered), function(i) {
    t <- exposure_prop_filtered$time[i]
    sids <- assign_sample_id(t)
    do.call(rbind, lapply(sids, function(sid) {
      data.frame(
        Sample_ID = sid,
        signature = exposure_prop_filtered$signature[i],
        exposure = exposure_prop_filtered$exposure[i],
        stringsAsFactors = FALSE
      )
    }))
  }))

  # Aggregate and reshape
  avg_exposure <- expanded_rows %>%
    dplyr::group_by(Sample_ID, signature) %>%
    dplyr::summarise(mean_exposure = mean(exposure), .groups = "drop")

  exposure_matrix <- avg_exposure %>%
    tidyr::pivot_wider(
      names_from = signature,
      values_from = mean_exposure,
      values_fill = list(mean_exposure = 0)
    ) %>%
    dplyr::arrange(Sample_ID)

  # Keep only SBS signatures
  exposure_sbs <- exposure_matrix %>%
    dplyr::select(Sample_ID, starts_with("SBS"))

  return(exposure_sbs)
}


load_signature_data <- function(paths) {
  lapply(paths, function(p) {
    if (!file.exists(p)) {
      warning("File not found: ", p)
      return(NULL)
    }

    if (grepl("\\.rds$", p)) {
      return(readRDS(p))
    } else if (grepl("\\.txt$", p)) {
      return(read.delim(p, header = TRUE, sep = "\t"))
    } else {
      warning(paste("Unknown file type:", p))
      return(NULL)
    }
  })
}



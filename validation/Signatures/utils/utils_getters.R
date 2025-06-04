qc_args <- function(name, coverage, purity, timepoint, sample_id, type) {
  #if (!is.numeric(name) || name %% 1 != 0) {
    #stop("Error: 'name' must be integer") # name : integer
  #}
  if (!is.numeric(coverage) || coverage %% 1 != 0) {
    stop("Error: 'coverage' must be integer") # coverage : integer
  }
  if (!is.numeric(purity) || (purity > 1 || purity < 0)) {
    stop("Error: 'purity' must be a number between 0 and 1") # purity : float [0, 1]
  }
  if (!is.null(timepoint) && (!is.numeric(timepoint) || timepoint %% 1 != 0)) {
    stop("Error: 'timepoint' must be integer") # timepoint : integer
  }
  if (!is.null(sample_id) && (!is.numeric(sample_id) || sample_id %% 1 != 0)) {
    stop("Error: 'sample_id' must be integer") # sample_id : integer
  }
  if (!(type %in% c("snv", "cna", "phylo")))    stop("Error: 'type' must be one of 'SNV', 'CNA' and 'phylo'")
 
  return(TRUE)
}



gen_ptr <- function(MAIN_PATH, SPN_ID, coverage, purity, timepoint, sample_id, type) {

  if (type == "cna") {
    path <- paste0(MAIN_PATH, SPN_ID, "/process/cna_data/")
    timepoint_ph <- ifelse(test = is.null(timepoint), yes = "\\d+", no = timepoint)
    sample_id_ph <- ifelse(test = is.null(sample_id), yes = "\\d+", no = sample_id)
    ptr <- paste0(SPN_ID, "_", timepoint_ph, "\\.", sample_id_ph, "_cna\\.rds$")
    return( c(path, ptr) )
  }

  else if (type == "snv") {
    path <- paste0(MAIN_PATH, SPN_ID, "/process/purity_", purity, "/")
    coverage_ph <- ifelse(test = is.null(coverage), yes = "\\d+", no = coverage)
    #purity_ph <- ifelse(test = is.null(purity), yes = "\\d+", no = purity)
    ptr <- paste0("seq_results_muts_merged_coverage_", coverage_ph, "x.rds")
    return( c(path, ptr) )
  }

  else if (type == "phylo") {
    path <- paste0(MAIN_PATH, SPN_ID, "/process/")
    ptr <- paste0("phylo_forest.sff")
    return( c(path, ptr) )
  }

  else {
    stop("Error in 'gen_ptr' function")
  }
}




getter_process <- function(MAIN_PATH = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/", 
                         SPN_ID = "SPN03", 
                         coverage = 100, 
                         purity = 0.6, 
                         timepoint = NULL, 
                         sample_id = NULL, 
                         type = "phylo") {
  type <- tolower(type)
  
  if ( !endsWith(MAIN_PATH, "/") ) {
    MAIN_PATH <- paste0(MAIN_PATH, "/")
  }
  
  tryCatch(
    {
      qc_args(SPN_ID, coverage, purity, timepoint, sample_id, type)
      
      #files_list <- search_files(directory=PATH, ptr)
      
      if ( type %in% c("snv", "cna") ) {
        x <- gen_ptr(MAIN_PATH, SPN_ID, coverage, purity, timepoint, sample_id, type)
        path <- x[1]
        ptr <- x[2]
        files_list <- search_files(path, ptr)
        if (length(files_list) != 0) {
          return(load_multiple_rds(files_list))
        }
      }
      else if (type=="phylo") {
        x <- gen_ptr(MAIN_PATH, SPN_ID, coverage, purity, timepoint, sample_id, type)
        #f <- paste0(PATH, SPN_ID, "/process/phylo_forest.sff")
        return(ProCESS::load_phylogenetic_forest(paste0(x[1], x[2])))
      }
      
      # if nothing found
      else {
        message("NOT FOUND!")
        return(NULL)
      }
      
    },
    error = function(cond) {
      message(conditionMessage(cond))
    }, 
    warning = function(cond) {
      message(conditionMessage(cond))
    }, 
    finally = {
    }
  )
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


get_sbs_exposures <- function(phylo_forest, sample_forest_path, snapshot_path) {

  # Load data
  exposure_prop <- phylo_forest$get_exposures()
  samples_forest <- load_samples_forest(sample_forest_path)
  samples_data <- samples_table(snapshot = snapshot_path, sample_forest = samples_forest)

  sample_ids <- samples_data %>%
    dplyr::filter(Samping_Time == 0) %>%
    dplyr::pull(Sample_ID)

  # Conditionally filter out time = 0 exposures
  if (any(exposure_prop$time != 0)) {
    exposure_prop_filtered <- exposure_prop %>% dplyr::filter(time != 0)
  } else {
    exposure_prop_filtered <- exposure_prop
  }

  # Function to assign sample ID based on time
  assign_sample_id <- function(t) {
    if (t == 0) {
      return(sample_ids)
    } else {
      match_idx <- which.min(abs(samples_data$Samping_Time - t))
      return(samples_data$Sample_ID[match_idx])
    }
  }

  # Expand timepoint exposures to sample IDs
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

  # Reshape
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


get_sigprofiler <- function(dataset, context) {

  if (context == "SBS") {
    ctx = "SBS96"
  } else if (context == "DBS") {
    ctx = "DBS78"
  } else {
    stop("Wrong context!")
  }

  base_path <- paste0("signature_deconvolution/SigProfiler/", dataset, "/results/", ctx, "/Suggested_Solution/")

  output <- list(
    COSMIC_EXPOSURE_PATH=paste0(base_path, "COSMIC_", ctx, "_Decomposed_Solution/Activities/COSMIC_", ctx, "_Activities.txt"),
    COSMIC_SIGNATURES_PATH=paste0(base_path, "COSMIC_", ctx, "_Decomposed_Solution/Signatures/COSMIC_", ctx, "_Signatures.txt"),
    DENOVO_EXPOSURE_PATH=paste0(base_path, ctx, "_De-Novo_Solution/Activities/", ctx, "_De-Novo_Activities_refit.txt"),
    DENOVO_SIGNATURE_PATH=paste0(base_path, ctx, "_De-Novo_Solution/Signatures/", ctx, "_De-Novo_Signatures.txt"),
    mut_counts=paste0("signature_deconvolution/SigProfiler/", dataset, "/results/", ctx, "Samples.txt")
  )
  return(output)
}



get_sparsesignatures <- function(dataset) {

  base_path <- paste0("signature_deconvolution/SparseSignatures/", dataset, "/")

  output <- list(
    NMF=paste0(base_path, dataset, "_nmf_Lasso_out.rds")
  )
  return(output)
}



get_signatures <- function(PATH, SPN, coverage, purity, dataset, context) {
  main_path <- paste0(PATH, SPN, "/tumourevo/", coverage, "x_", purity, "p_mutect2_ascat/")

  sparse_paths <- get_sparsesignatures(dataset)
  sigprofiler_paths <- get_sigprofiler(dataset, context)

  full_paths <- c(
    lapply(sparse_paths, function(p) paste0(main_path, p)),
    lapply(sigprofiler_paths, function(p) paste0(main_path, p))
  )

  return(unlist(full_paths))
}



load_signature_data <- function(paths) {
  lapply(paths, function(p) {
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

setwd("/orfeo/cephfs/scratch/cdslab/kdavydzenka/ProCESS-examples/validation/Signatures")

pkgs <- c("tidyverse", "ggplot2", "caret", "ggtext", "reshape2", "lsa", "Metrics", "MutationalPatterns",
"ggalluvial")

sapply(pkgs, require, character.only = TRUE)

source("utils/utils_getters.R")
source("utils/utils_sparsesig.R")
source("utils/utils_validation.R")
source("utils/utils.R")

### Get ProCESS exposure data ###

main_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
spn_list <- c("SPN01", "SPN02", "SPN03", "SPN04")
coverage_list <- c(50)
purity_list <- c(0.3, 0.6)

process_exposure <- list()

for (spn in spn_list) {
  for (cov in coverage_list) {
    for (pur in purity_list) {

      # Get phylo_forest
      phylo_forest <- getter_process(
        MAIN_PATH = main_path,
        SPN_ID = spn,
        coverage = cov,
        purity = pur,
        timepoint = NULL,
        sample_id = NULL,
        type = "phylo"
      )

      # Build dynamic file paths
      sample_forest_path <- file.path(spn, "process", "sample_forest.sff")
      snapshot_path <- file.path(spn, "process", spn)

      # Try to extract exposures
      exposure_sbs <- tryCatch({
        get_sbs_exposures(
          phylo_forest,
          sample_forest_path = sample_forest_path,
          snapshot_path = snapshot_path
        )
      }, error = function(e) {
        warning("Failed for ", spn, " at cov=", cov, ", pur=", pur, ": ", e$message)
        return(NULL)
      })

      # Save to results
      results_list[[paste(spn, cov, pur, sep = "_")]] <- exposure_sbs
    }
  }
}


### Get tumourevo exposure data ###

base_path <- ""
dataset <- "SCOUT"
context <- "SBS"

tumourevo_signature_res <- list()

for (spn in spn_list) {
  tumourevo_signature_res[[spn]] <- list()

  for (cov in coverage_list) {
    tumourevo_signature_res[[spn]][[paste0("coverage_", cov)]] <- list()

    for (pur in purity_list) {
      message("Processing ", spn, " | Coverage: ", cov, " | Purity: ", pur)

      result <- tryCatch({
        paths <- get_signatures(
          PATH = base_path,
          SPN = spn,
          coverage = cov,
          purity = pur,
          dataset = dataset,
          context = context
        )

        data <- load_signature_data(paths)

        names(data) <- c(
          "Sparsesig",
          "Sigprofiler_CosmicExposure",
          "Sigprofiler_CosmicSignature",
          "Sigprofiler_DenovoExposure",
          "Sigprofiler_DenovoSignature",
          "mut_counts"
        )

        data
      }, error = function(e) {
        message("Failed for ", spn, " | cov: ", cov, " | pur: ", pur)
        message("Reason: ", e$message)
        NULL
      })

      # Store only if result was successful
      if (!is.null(result)) {
        tumourevo_signature_res[[spn]][[paste0(cov)]][[paste0("p", pur)]] <- result
      }
    }
  }
}


### Map De novo SparseSignatures results to COSMIC using cosine similarity ###

cosmic_path <- "COSMIC_v3.4/COSMIC_v3.4_SBS_GRCh38.txt"
sparsig_cosmic <- list()

for (spn in names(tumourevo_signature_res)) {

  # Set threshold depending on SPN
  threshold <- switch(spn,
                      "SPN01" = 0.5,
                      "SPN03" = 0.7,
                      "SPN04" = 0.7,
                      0.7)  # default if any others

  for (coverage in names(tumourevo_signature_res[[spn]])) {
    for (purity in names(tumourevo_signature_res[[spn]][[coverage]])) {

      # Access Sparsesig output and mutation counts
      sparsesig_out <- tumourevo_signature_res[[spn]][[coverage]][[purity]][["Sparsesig"]]
      mut_counts <- tumourevo_signature_res[[spn]][[coverage]][[purity]][["mut_counts"]]

      # Map to COSMIC signatures with appropriate threshold
      remapped_exposures_prop <- map_sparsesig_to_cosmic(
        sparsesig_out = sparsesig_out,
        mut_counts = mut_counts,
        cosmic_path = cosmic_path,
        threshold = threshold
      )

      # Store results
      sparsig_cosmic[[spn]][[coverage]][[purity]] <- remapped_exposures_prop
    }
  }
}


### Validate Signatures across combinations ###

spn_list <- c("SPN01", "SPN03", "SPN04")
coverage <- 50
purity <- 0.6

ground_truth <- list()

for (spn in spn_list) {
  if (spn %in% names(process_exposures_list)) {
    gt <- process_exposures_list[[spn]] %>%
      tibble::column_to_rownames("Sample_ID") %>%
      as.matrix()
    
    ground_truth[[spn]] <- gt
  } else {
    warning(paste("Missing ground truth for", spn))
  }
}

sparsesig_aligned <- align_sparsesig_res(sparsig_cosmic)
sigprof_aligned <- align_sigprofiler_res(tumourevo_signature_res)


ground_truth_nested <- list()

for (spn in names(ground_truth)) {
  ground_truth_nested[[spn]] <- list()
  coverage_key <- paste0("coverage_", coverage)  # "coverage_50"
  purity_key <- paste0("purity_", purity)        # "purity_0.6"
  
  ground_truth_nested[[spn]][[coverage_key]] <- list()
  ground_truth_nested[[spn]][[coverage_key]][[purity_key]] <- ground_truth[[spn]]
}

metrics_sparsesig <- evaluate_all_combined(ground_truth_nested, sparsesig_aligned, threshold = 0.05)
metrics_sigprof <- evaluate_all_combined(ground_truth_nested, sigprof_aligned, threshold = 0.05)

combined_metrics <- bind_rows(
  metrics_sparsesig %>% mutate(Tool = "SparseSignatures"),
  metrics_sigprof %>% mutate(Tool = "SigProfiler")
)

saveRDS(combined_metrics, file = "combined_metrics_signatures.rds")


# Sankey plot - Compare estimated and true signatures  #

sankey_df <- prepare_sankey_data(ground_truth_nested, sparsesig_aligned, sigprof_aligned)

sankey_df <- sankey_df %>%
  dplyr::mutate(
    Coverage = as.numeric(gsub("coverage_", "", Coverage)),
    Purity = as.numeric(gsub("purity_", "", Purity))
  )

spns <- c("SPN01", "SPN03", "SPN04")
sankey_plots <- lapply(spns, function(spn) generate_sankey(spn, cov = 50, pur = 0.6))

wrapped_sankey <- wrap_plots(sankey_plots, ncol = 3)
wrapped_sankey


### Exposure validation ###

setwd("/orfeo/cephfs/scratch/cdslab/kdavydzenka/ProCESS-examples/validation/Signatures")

pkgs <- c("tidyverse", "ggplot2", "caret", "ggtext", "reshape2", "lsa", "Metrics", "pheatmap", "MutationalPatterns", "ggvenn")
sapply(pkgs, require, character.only = TRUE)

source("utils/utils_getters.R")
source("utils/utils_validation.R")

### Get ProCESS exposure data ###

main_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
spn_list <- c("SPN01", "SPN02", "SPN03", "SPN04")
coverage_list <- c(50)
purity_list <- c(0.3, 0.6)

results_list <- list()

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



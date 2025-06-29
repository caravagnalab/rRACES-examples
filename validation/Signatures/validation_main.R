setwd("/orfeo/cephfs/scratch/cdslab/kdavydzenka/ProCESS-examples/")

pkgs <- c("tidyverse", "ggplot2", "caret", "ggtext", "reshape2", "lsa", "Metrics", "MutationalPatterns",
"ggalluvial", "patchwork")

sapply(pkgs, require, character.only = TRUE)

source("getters/process_getters.R")
source("getters/tumourevo_getters.R")
source("validation/signatures/utils/utils_getters.R")
source("validation/signatures/utils/utils_sparsesig.R")
source("validation/signatures/utils/utils_validation.R")
source("validation/signatures/utils/utils.R")
source("validation/signatures/utils/utils_plots.R")

### Get ProCESS exposure data ###

base_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
spn_id <- c("SPN01", "SPN03", "SPN04")

process_exposures_list <- list()

for (spn_id in spn) {
  exposure_sbs <- tryCatch({
    get_sbs_exposures(spn = spn_id, base_path = base_path)
  }, error = function(e) {
    warning("Failed for ", spn_id, ": ", e$message)
    return(NULL)
  })

  process_exposures[[spn_id]] <- exposure_sbs
}


### Get tumourevo data ###

spn_list <- c("SPN01", "SPN03", "SPN04")
coverage_list <- c(50)
purity_list <- c(0.3, 0.6, 0.9)
dataset <- "SCOUT"
context <- "SBS96"
vcf_caller <- "mutect2"
cna_caller <- "ascat"


tumourevo_signature_res <- list()

# Iterate over all combinations
for (spn in spn_list) {
  tumourevo_signature_res[[spn]] <- list()
  
  for (cov in coverage_list) {
    tumourevo_signature_res[[spn]][[paste0("coverage_", cov)]] <- list()
    
    for (pur in purity_list) {
      message("Processing ", spn, " | Coverage: ", cov, " | Purity: ", pur)
      
      result <- tryCatch({
        # Get SparseSignatures paths
        sparsesig <- get_tumourevo_signatures(
          spn = spn,
          coverage = cov,
          purity = pur,
          vcf_caller = vcf_caller,
          cna_caller = cna_caller,
          tool = "SparseSignatures",
          base_path = base_path
        )
        
        # Get SigProfiler paths
        sigprofiler <- get_tumourevo_signatures(
          spn = spn,
          coverage = cov,
          purity = pur,
          vcf_caller = vcf_caller,
          cna_caller = cna_caller,
          tool = "SigProfiler",
          context = context,  # SBS96!
          base_path = base_path
        )
        
        # Combine and load paths
        paths <- c(
          sparsesig$nmf_Lasso_out,
          sparsesig$cv_means_mse,
          sparsesig$best_params_config,
          sigprofiler$COSMIC_exposure,
          sigprofiler$COSMIC_signatures,
          sigprofiler$denovo_exposure,
          sigprofiler$denovo_signatures
        )
        
        data <- load_signature_data(paths)
        
        names(data) <- c(
          "SparseSig_nmf_Lasso_out",
          "SparseSig_cv_means_mse",
          "SparseSig_best_params_config",
          "SigProfiler_COSMIC_exposure",
          "SigProfiler_COSMIC_signatures",
          "SigProfiler_denovo_exposure",
          "SigProfiler_denovo_signatures"
        )
        
        data
      }, error = function(e) {
        message("Failed for ", spn, " | cov: ", cov, " | pur: ", pur)
        message("Reason: ", e$message)
        NULL
      })
      
      # Store result if successful
      if (!is.null(result)) {
        tumourevo_signature_res[[spn]][[paste0("coverage_", cov)]][[paste0("purity_", pur)]] <- result
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
      sparsesig_out <- tumourevo_signature_res[[spn]][[coverage]][[purity]][["SparseSig_nmf_Lasso_out"]]
      mut_counts <- tumourevo_signature_res[[spn]][[coverage]][[purity]][["SigProfiler_mut_counts"]]

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
purity_values <- c(0.3, 0.6, 0.9)

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
  coverage_key <- paste0("coverage_", coverage)
  
  if (is.null(ground_truth_nested[[spn]])) {
    ground_truth_nested[[spn]] <- list()
  }
  
  if (is.null(ground_truth_nested[[spn]][[coverage_key]])) {
    ground_truth_nested[[spn]][[coverage_key]] <- list()
  }
  
  for (purity in purity_values) {
    purity_key <- paste0("purity_", purity)
    ground_truth_nested[[spn]][[coverage_key]][[purity_key]] <- ground_truth[[spn]]
  }
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
  ) %>% 
  dplyr::filter(Coverage == 50, Purity == 0.6)



spns <- c("SPN01", "SPN03", "SPN04")
sankey_plots <- lapply(spns, function(spn) generate_sankey(spn, cov = 50, pur = 0.6))
wrapped_sankey <- wrap_plots(sankey_plots, ncol = 3)
wrapped_sankey


## Exposure validation ##

# Align exposure data

aligned_exposures <- list(
  SparseSignatures = list(),
  SigProfiler = list()
)

spn_list <- names(ground_truth_nested)

for (spn in spn_list) {
  coverage_list <- names(ground_truth_nested[[spn]])

  for (coverage in coverage_list) {
    purity_list <- names(ground_truth_nested[[spn]][[coverage]])

    for (purity in purity_list) {
      # Access exposure matrices
      gt <- ground_truth_nested[[spn]][[coverage]][[purity]]
      sparse <- sparsesig_aligned[[spn]][[coverage]][[purity]]
      sigprof <- sigprof_aligned[[spn]][[coverage]][[purity]]

      # Find shared signatures
      shared_sparse <- intersect(colnames(gt), colnames(sparse))
      shared_sigprof <- intersect(colnames(gt), colnames(sigprof))

      # Subset to shared signatures
      gt_sparse <- gt[, shared_sparse, drop = FALSE]
      sparse <- sparse[, shared_sparse, drop = FALSE]

      gt_sigprof <- gt[, shared_sigprof, drop = FALSE]
      sigprof <- sigprof[, shared_sigprof, drop = FALSE]

      # Create unique key for storage
      key <- paste(spn, coverage, purity, sep = "_")

      # Store aligned pairs
      aligned_exposures$SparseSignatures[[key]] <- list(gt = gt_sparse, tool = sparse)
      aligned_exposures$SigProfiler[[key]]     <- list(gt = gt_sigprof, tool = sigprof)
    }
  }
}


# Calculate Cosine similarity and MSE

compute_exposure_metrics <- function(aligned_exposures) {
  cosine_results <- list()
  mse_results <- list()
  cosine_idx <- 1
  mse_idx <- 1

  for (method in names(aligned_exposures)) {
    keys <- names(aligned_exposures[[method]])

    for (key in keys) {
      # Extract matrices
      gt_mat <- aligned_exposures[[method]][[key]]$gt
      tool_mat <- aligned_exposures[[method]][[key]]$tool

      # Skip if ground truth is empty
      if (nrow(gt_mat) == 0 || ncol(gt_mat) == 0) next

      # Parse metadata from key: expected format "SPN01_coverage_50_purity_0.9"
      parts <- str_split(key, "_", simplify = TRUE)
      SPN <- parts[1]
      Coverage <- as.numeric(parts[3])
      Purity <- as.numeric(parts[5])

      ## --- COSINE SIMILARITY ---
      # Align sample rows
      common_samples <- intersect(rownames(gt_mat), rownames(tool_mat))
      if (length(common_samples) == 0) next

      gt_common <- as.matrix(gt_mat[common_samples, , drop = FALSE])
      tool_common <- as.matrix(tool_mat[common_samples, , drop = FALSE])

      # Compute cosine similarity
      sims <- cosine_similarity_exposures(gt_common, tool_common)

      cosine_results[[cosine_idx]] <- tibble(
        Sample = common_samples,
        SPN = SPN,
        Coverage = Coverage,
        Purity = Purity,
        Tool = method,
        CosineSimilarity = sims
      )
      cosine_idx <- cosine_idx + 1

      # Align signatures
      common_sigs <- intersect(colnames(gt_mat), colnames(tool_mat))
      if (length(common_sigs) == 0) next

      gt_sub <- as.matrix(gt_mat[, common_sigs, drop = FALSE])
      tool_sub <- as.matrix(tool_mat[, common_sigs, drop = FALSE])

      # Compute MSE
      mse_vec <- colMeans((gt_sub - tool_sub)^2, na.rm = TRUE)

      mse_results[[mse_idx]] <- tibble(
        SPN = SPN,
        Coverage = Coverage,
        Purity = Purity,
        Tool = method,
        Signature = names(mse_vec),
        MSE = mse_vec
      )
      mse_idx <- mse_idx + 1
    }
  }

  # Combine results into data frames
  list(
    cosine = bind_rows(cosine_results),
    mse = bind_rows(mse_results)
  )
}

metrics <- compute_exposure_metrics(aligned_exposures)



### Generate final report plot ###
p_cosine <- plot_cosine_similarity(metrics[["cosine"]])
p_mse <- plot_mse_per_signature(metrics[["mse"]])

bottom_plots <- p_cosine + p_mse +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

bottom_plots

final_plots <- wrapped_sankey / bottom_plots +
  plot_layout(heights = c(1, 1))

final_plots

ggsave("final_plot.pdf", plot = final_plots,
       width = 20, height = 10, units = "in")




   

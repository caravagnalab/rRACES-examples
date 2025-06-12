# Reshape data

reshape_exposures_long <- function(exposures_mat, spn, coverage = NA, purity = NA, method_name) {
  df <- as.data.frame(exposures_mat)
  df$Sample_ID <- rownames(exposures_mat)
  
  df_long <- df %>%
    pivot_longer(-Sample_ID, names_to = "Signature", values_to = "Exposure") %>%
    mutate(SPN = spn,
           Coverage = coverage,
           Purity = purity,
           Method = method_name)
  
  return(df_long)
}

extract_ground_truth_long <- function(ground_truth_list) {
  out <- list()

  for (spn in names(ground_truth_list)) {
    for (coverage in names(ground_truth_list[[spn]])) {
      for (purity in names(ground_truth_list[[spn]][[coverage]])) {
        exposures_mat <- ground_truth_list[[spn]][[coverage]][[purity]]

        if (is.null(exposures_mat)) next

        long_df <- reshape_exposures_long(exposures_mat, spn, coverage, purity, "GroundTruth")
        out[[paste(spn, coverage, purity, sep = "_")]] <- long_df
      }
    }
  }

  do.call(rbind, out)
}


extract_sparsesig_long <- function(sparsesig_aligned) {
  out <- list()

  for (spn in names(sparsesig_aligned)) {
    for (coverage in names(sparsesig_aligned[[spn]])) {
      for (purity in names(sparsesig_aligned[[spn]][[coverage]])) {
        df <- sparsesig_aligned[[spn]][[coverage]][[purity]]

        if (is.null(df)) next

        # df is samples x signatures matrix/data.frame, samples are rownames
        long_df <- reshape_exposures_long(df, spn, coverage, purity, "SparseSignatures")
        out[[paste(spn, coverage, purity, sep = "_")]] <- long_df
      }
    }
  }

  do.call(rbind, out)
}


extract_sigprofiler_long <- function(sigprof_aligned) {
  out <- list()

  for (spn in names(sigprof_aligned)) {
    for (coverage in names(sigprof_aligned[[spn]])) {
      for (purity in names(sigprof_aligned[[spn]][[coverage]])) {
        df <- sigprof_aligned[[spn]][[coverage]][[purity]]

        if (is.null(df)) next

        long_df <- reshape_exposures_long(df, spn, coverage, purity, "SigProfiler")
        out[[paste(spn, coverage, purity, sep = "_")]] <- long_df
      }
    }
  }

  do.call(rbind, out)
}

prepare_sankey_data <- function(ground_truth_list, sparsesig_aligned, sigprof_aligned) {
  gt_long <- extract_ground_truth_long(ground_truth_list)
  sparsesig_long <- extract_sparsesig_long(sparsesig_aligned)
  sigprof_long <- extract_sigprofiler_long(sigprof_aligned)

  combined <- bind_rows(gt_long, sparsesig_long, sigprof_long)

  # Filter to keep only exposures > 0
  combined <- combined %>% filter(Exposure > 0)

  return(combined)
}




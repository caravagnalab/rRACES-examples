require(patchwork)
require(RColorBrewer)
require(ggExtra)
require(ggVennDiagram)
require(ggalluvial)
require(ggplotify)

string_difference <- function(str1, str2) {
  str1_chars <- unlist(strsplit(str1, ""))
  str2_chars <- unlist(strsplit(str2, ""))
  
  for (char in str2_chars) {
    match_idx <- match(char, str1_chars)
    if (!is.na(match_idx)) {
      str1_chars <- str1_chars[-match_idx]
    }
  }
  
  return(paste(str1_chars, collapse = ""))
}

merge_datasets <- function(snp_caller, ground_truth, min_vaf) {
  df <- snp_caller %>%
    full_join(ground_truth, by = c("mutationID"), suffix = c("_caller", "_truth")) %>%
      dplyr::mutate(
        positive_truth = !is.na(VAF_truth) & VAF_truth > min_vaf,
        positive_call = !is.na(VAF_caller) & VAF_caller > min_vaf,
        true_positive = positive_truth & positive_call,
        false_positive = !positive_truth & positive_call,
        false_negative = positive_truth & !positive_call,
        true_negative = !positive_truth & !positive_call
      ) %>% 
      dplyr::mutate(VAF_truth = ifelse(is.na(VAF_truth), 0, VAF_truth)) %>% 
      dplyr::mutate(VAF_caller = ifelse(is.na(VAF_caller), 0, VAF_caller)) %>% 
      dplyr::mutate(DP_truth = ifelse(is.na(DP_truth), 0, DP_truth)) %>% 
      dplyr::mutate(DP_caller = ifelse(is.na(DP_caller), 0, DP_caller))
  
  return(df)
}

get_colors = function(df) {
  all_levels <- levels(factor(df$FILTER))
  full_colors <- setNames(RColorBrewer::brewer.pal(8, "Set2")[1:length(all_levels)], all_levels)  
  
  original_names = names(full_colors)
  if ("PASS" %in% names(full_colors)) {
    names(full_colors)[which(names(full_colors)=="PASS")] = original_names[1]
    names(full_colors)[1] = "PASS"
  }
  
  if ("Other" %in% names(full_colors)) {
    names(full_colors)[which(names(full_colors)=="Other")] = original_names[2]
    names(full_colors)[2] = "Other"
  }
  
  full_colors
}

get_sample_info <- function(vcf_path) {
  # Extract caller name from vcf_path
  caller_match <- regexpr("strelka|races|freebayes|mutect2|haplotypecaller", vcf_path)
  caller_name <- ifelse(caller_match != -1, regmatches(vcf_path, caller_match), NA)
  
  # Extract sample name from vcf_path
  sample_match <- regmatches(vcf_path, regexpr("SPN[0-9]+_[0-9]+\\.[0-9]+", vcf_path))
  sample_name <- ifelse(length(sample_match) > 0, sample_match, NA)
  
  list(sample_name = sample_name, caller_name = caller_name)
}

get_report <- function(seq_res_long, caller_res, sample_info, min_vaf) {
  # Group filter categories for better visualization
  caller_res <- caller_res %>% 
    dplyr::group_by(FILTER) %>% 
    dplyr::mutate(f = n() / nrow(caller_res)) %>% 
    dplyr::mutate(FILTER = ifelse(FILTER == "PASS" | f > .05, FILTER, "Other"))
  
  # Create merged dataset for all mutations
  merged_df <- merge_datasets(caller_res, seq_res_long, min_vaf)
  
  # Keep only PASS mutations
  caller_res_pass <- caller_res %>% 
    dplyr::filter(FILTER == "PASS" & !is.na(FILTER))
  
  # Create merged dataset with 0 threshold for PASS mutations
  merged_df0 <- merge_datasets(caller_res_pass, seq_res_long, 0)
  
  # Plot false negative VAF distribution
  p_false_negative_VAF_dist <- plot_false_negative_vaf_dist(merged_df0)
  
  # Plot basic histograms
  races_coverage <- plot_races_coverage(seq_res_long)
  races_VAF <- plot_races_vaf(seq_res_long)
  caller_coverage <- plot_caller_coverage(caller_res_pass, sample_info)
  caller_VAF <- plot_caller_vaf(caller_res_pass, sample_info)
  
  # Plot VAF and coverage differences
  vaf_differences <- plot_vaf_difference(merged_df)
  cov_differences <- plot_cov_difference(merged_df)
  
  # Get color palette for FILTER categories
  colors <- get_colors(merged_df)
  
  # Plot filter distribution
  p_filter_dist <- plot_filter_distribution(merged_df, colors)
  
  # --- All mutations ---
  # Plot scatter plots for all mutations
  p_scatter_VAF_all <- plot_vaf_scatter_all(merged_df, colors, sample_info)
  p_scatter_DP_all <- plot_dp_scatter_all(merged_df, colors, sample_info)
  
  # Calculate performance metrics for all mutations
  y_true <- as.numeric(factor((merged_df$VAF_truth > min_vaf), levels=c(FALSE, TRUE))) - 1
  y_pred <- as.numeric(factor((merged_df$VAF_caller > min_vaf), levels=c(FALSE, TRUE))) - 1
  
  # Plot Venn diagram for all mutations
  p_venn_all <- plot_venn_diagram(merged_df, sample_info$caller_name) +
    ggtitle("All called mutations")
  
  # Calculate metrics for all mutations
  metrics_all <- compute_metrics_from_vectors(y_true, y_pred)
  
  # --- PASS mutations ---
  # Create merged dataset for PASS mutations
  merged_df_pass <- merge_datasets(caller_res_pass, seq_res_long, min_vaf)
  
  # Plot scatter plots for PASS mutations
  p_scatter_VAF_pass <- plot_vaf_scatter_pass(merged_df_pass, colors, sample_info)
  p_scatter_DP_pass <- plot_dp_scatter_pass(merged_df_pass, colors, sample_info)
  
  # Calculate performance metrics for PASS mutations
  y_true_pass <- as.numeric(factor((merged_df_pass$VAF_truth > min_vaf), levels=c(FALSE, TRUE))) - 1
  y_pred_pass <- as.numeric(factor((merged_df_pass$VAF_caller > min_vaf), levels=c(FALSE, TRUE))) - 1
  
  # Plot Venn diagram for PASS mutations
  p_venn_pass <- plot_venn_diagram(merged_df_pass, sample_info$caller_name) +
    ggtitle("Only PASS mutations")
  
  # Calculate metrics for PASS mutations
  metrics_pass <- compute_metrics_from_vectors(y_true_pass, y_pred_pass)
  
  # Combine metrics for both sets
  metrics <- dplyr::bind_rows(
    metrics_all %>% tidyr::pivot_longer(cols = colnames(metrics_all)) %>% 
      dplyr::mutate(Mutations = "All"),
    metrics_pass %>% tidyr::pivot_longer(cols = colnames(metrics_pass)) %>% 
      dplyr::mutate(Mutations = "Only Pass")  
  )
  
  # Plot metrics comparison
  p_metrics <- plot_metrics_comparison(metrics)
  
  metrics_over_VAF = plot_metric_over_VAF_threshold(seq_res_long, caller_res, only_pass = TRUE)
  
  # VAF comparison for true positives
  vaf_comparison <- merged_df_pass %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(true_positive) %>% 
    dplyr::group_by(chr_caller) %>% 
    dplyr::summarise(cor_coeff = cor.test(VAF_caller, VAF_truth)$estimate, RMSE = RMSE(VAF_caller, VAF_truth))
  
  # Design layout for combined report
  design <- "
  AABBCC
  AABBCC
  DDEEFF
  DDEEFF
  GGHHII
  GGHHII
  GGHHII
  LLMMNN
  LLMMNN
  LLMMNN
  OOOPPP
  OOOPPP
  "
  
  # Generate title and subtitle for the report
  title <- paste0(sample_info$spn, ", sample ", sample_info$sample_id, 
                  ", calls by ", sample_info$caller_name)
  subtitle <- paste0("Only ", sample_info$mut_type, " mutations, purity = ", 
                     sample_info$purity, ", coverage = ", sample_info$coverage, "x")
  
  # Combine all plots into a single report
  report_plot <- free(races_coverage) + free(caller_coverage) + free(p_filter_dist) +
    free(races_VAF) + free(caller_VAF) + free(p_false_negative_VAF_dist) +
    free(p_scatter_DP_all) + free(p_scatter_VAF_all) + free(p_venn_all) +
    free(p_scatter_DP_pass) + free(p_scatter_VAF_pass) + free(p_venn_pass) +
    free(p_metrics) + free(metrics_over_VAF) +
    plot_layout(design = design) +
    plot_annotation(title, subtitle) & 
    theme(text = element_text(size = 12))
  
  # Return results
  list(
    report_plot = report_plot, 
    report_metrics = metrics, 
    vaf_comparison = vaf_comparison
  )
}


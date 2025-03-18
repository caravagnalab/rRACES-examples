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

merge_datasets <- function(snp_caller, ground_truth) {
  df <- snp_caller %>%
    full_join(ground_truth, by = c("mutationID"), suffix = c(".caller", ".races")) %>%
    mutate(
      positive_truth = !is.na(BAF.races) & BAF.races > 0,
      positive_call = !is.na(BAF.caller) & BAF.caller > 0,
      true_positive = positive_truth & positive_call,
      false_positive = !positive_truth & positive_call,
      false_negative = positive_truth & !positive_call,
      true_negative = !positive_truth & !positive_call
    ) %>% 
    dplyr::mutate(BAF.races = ifelse(is.na(BAF.races), 0, BAF.races)) %>% 
    dplyr::mutate(BAF.caller = ifelse(is.na(BAF.caller), 0, BAF.caller)) %>% 
    dplyr::mutate(DP.races = ifelse(is.na(DP.races), 0, DP.races)) %>% 
    dplyr::mutate(DP.caller = ifelse(is.na(DP.caller), 0, DP.caller))
  
  return(df)
}


plot_venn_diagram = function(merged_df, caller_name) {
  id_muts_races = merged_df %>% dplyr::filter(positive_truth) %>% dplyr::pull(mutationID)
  id_muts_caller = merged_df %>% dplyr::filter(positive_call) %>% dplyr::pull(mutationID)
  x = list("races"=id_muts_races, caller_name=id_muts_caller)
  names(x) = c("races", caller_name)
  ggVennDiagram::ggVennDiagram(x) +
    scale_fill_gradient2(low = "white", high = "#4981BF", mid = "white", midpoint=0) +
    coord_flip() +
    theme(legend.position = "none")
}


seq_to_long <- function(seq_results) {
  # Extract sample names from column names
  sample_names <- strsplit(colnames(seq_results)[grepl(".VAF", colnames(seq_results), fixed = TRUE)], ".VAF") %>% unlist()
  
  seq_df <- lapply(sample_names, function(sn) {
    cc <- c("chr", "chr_pos", "ref", "alt", "causes", "classes", colnames(seq_results)[grepl(paste0(sn, "."), colnames(seq_results), fixed = TRUE)])
    seq_results[, cc] %>%
      `colnames<-`(c("chr", "chr_pos", "ref", "alt", "causes", "classes", "occurences", "coverage", "VAF")) %>%
      dplyr::mutate(sample_name = sn)
  }) %>% do.call("bind_rows", .)
  
  seq_df %>%
    dplyr::rename(chr = chr, from = chr_pos, DP = coverage, NV = occurences, ALT = alt) %>%
    dplyr::mutate(to = from)
}

plot_baf_difference <- function(df) {
  df <- df %>%
    mutate(BAF_diff = BAF.caller - BAF.races)
  
  ggplot(df, mapping = aes(x=FILTER, y=BAF_diff)) +
    geom_violin() +
    labs(title = "Boxplot of BAF differences", y = "BAF Difference (Caller - Truth)", x = "") +
    theme_bw()
}

# Function to plot Histogram of depth Differences
plot_cov_difference <- function(df) {
  df <- df %>%
    dplyr::mutate(DP_diff = DP.caller - DP.races)
  
  ggplot(df, mapping = aes(x=FILTER, y=DP_diff)) +
    geom_violin() +
    labs(title = "Boxplot of DP differences", y = "DP Difference (Caller - Truth)", x = "") +
    theme_bw()
}

plot_scatter_with_corr <- function(data, x_var, y_var, col_var = "FILTER", title = NULL) {
  # Perform correlation test to get r and p-value
  data[[x_var]][is.na(data[[x_var]])] = 0
  data[[y_var]][is.na(data[[y_var]])] = 0
  
  corr_test <- cor.test(data[[x_var]], data[[y_var]], use = "complete.obs")
  corr_coef <- corr_test$estimate
  p_val <- corr_test$p.value
  
  # Format p-value appropriately
  if (p_val < 0.001) {
    p_text <- "p < 0.001"
  } else if (p_val < 0.01) {
    p_text <- paste0("p = ", format(round(p_val, 3), nsmall = 3))
  } else if (p_val < 0.05) {
    p_text <- paste0("p = ", format(round(p_val, 2), nsmall = 2))
  } else {
    p_text <- paste0("p = ", format(round(p_val, 2), nsmall = 2))
  }
  
  # Format correlation coefficient (rounded to 2 decimal places)
  corr_text <- paste0("r = ", round(corr_coef, 2), ", ", p_text)
  
  # Create plot
  p <- ggplot(data, aes_string(x = x_var, y = y_var, col = col_var)) +
    geom_point(alpha = 0.5, size = .4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    annotate("text", x = Inf, y = -Inf, label = corr_text,
             hjust = 1.1, vjust = -0.5) +
    theme_bw() +
    labs(x = x_var, y = y_var)
  p
}


library(caret)

compute_metrics <- function(actual, predicted) {
  cm <- table(Actual = actual, Predicted = predicted)
  confusion_matrix <- confusionMatrix(as.factor(predicted), as.factor(actual), positive = "1")
  
  metrics <- dplyr::tibble(
    Accuracy = confusion_matrix$overall["Accuracy"],
    Sensitivity = confusion_matrix$byClass["Sensitivity"],
    Precision = confusion_matrix$byClass["Precision"],
    Recall = confusion_matrix$byClass["Recall"],
    F1_Score = confusion_matrix$byClass["F1"]
  )
  
  return(metrics)
}

plot_filter_distribution = function(merged_df, log_scale = TRUE) {
  p = merged_df %>% 
    dplyr::group_by(FILTER) %>% 
    dplyr::summarise(n=n()) %>% 
    na.omit() %>% 
    ggplot(mapping = aes(x=reorder(FILTER, +n), y=n)) +
    geom_col() +
    theme_bw() +
    coord_flip() +
    labs(x = "FILTER", y="count")  
  
  if (log_scale) {
    p <- p +
      scale_y_continuous(trans = "log10")
  }
  
  p
}



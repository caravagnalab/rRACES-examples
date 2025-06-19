library(caret)

method_colors = c(
  "ProCESS" = "gray80",
  "mutect2" = "lightsteelblue",
  "mutect2 (all)" = "steelblue4",  # darker steelblue
  "strelka" = "coral",
  "strelka (all)" = "coral4",     # darker coral
  "freebayes" = "#8FBC8B",
  "freebayes (all)" = "#228B22"   # darker green (ForestGreen)
)

# Function to compute confusion matrix and performance metrics
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

plot_venn_diagram = function(merged_df, caller_name) {
  id_muts_races = merged_df %>% dplyr::filter(positive_truth) %>% dplyr::pull(mutationID)
  id_muts_caller = merged_df %>% dplyr::filter(positive_call) %>% dplyr::pull(mutationID)
  x = list("races"=id_muts_races, caller_name=id_muts_caller)
  names(x) = c("races", caller_name)
  ggVennDiagram(x) +
    scale_fill_gradient2(low = "white", high = "#4981BF", mid = "white", midpoint=0) +
    coord_flip() +
    theme(legend.position = "none")
}

# Function to plot the confusion matrix
plot_confusion_matrix <- function(actual, predicted) {
  # Create the confusion matrix
  cm <- table(Actual = actual, Predicted = predicted)
  
  # Convert to data frame for ggplot
  cm_df <- as.data.frame(cm)
  
  # Calculate global percentages (out of total observations)
  total <- sum(cm)
  cm_df$Percentage <- (cm_df$Freq / total) * 100
  
  # Create labels with both count and percentage
  cm_df$Label <- paste0(cm_df$Freq, "\n(", round(cm_df$Percentage, 1), "%)")
  
  # Add a column to identify diagonal vs. off-diagonal elements
  cm_df$Diagonal <- ifelse(cm_df$Actual == cm_df$Predicted, "Correct", "Incorrect")
  
  # Plot
  ggplot(cm_df, aes(x = Predicted, y = Actual, fill = Diagonal)) +
    geom_tile() +
    geom_text(aes(label = Label), color = "black", size = 4) +
    scale_fill_manual(values = c("Correct" = alpha("#2E8B57", .9), "Incorrect" = alpha("indianred", .9))) +
    theme_bw() +
    labs(title = "Confusion Matrix", x = "Predicted", y = "Actual")
}

plot_scatter_with_corr <- function(data, x_var, y_var, col_var = "FILTER", title = NULL) {
  # Perform correlation test to get r and p-value
  data[[x_var]][is.na(data[[x_var]])] = 0
  data[[y_var]][is.na(data[[y_var]])] = 0
  #print(data[[x_var]])
  #print(data[[y_var]])
  #print(x_var)
  #print(y_var)
  if (nrow(data)<3){
    p <- ggplot(data, aes_string(x = x_var, y = y_var, col = col_var)) +
      geom_point(alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      theme_bw() +
      labs(x = x_var, y = y_var)
    return(p)
  }
  corr_test <- cor.test(data[[x_var]], data[[y_var]], use = "complete.obs")
  corr_coef <- corr_test$estimate
  p_val <- corr_test$p.value
  #print(corr_test)
  # Format p-value appropriately
  # if (p_val < 0.001) {
  #  p_text <- "p < 0.001"
  # } else if (p_val < 0.01) {
  #  p_text <- paste0("p = ", format(round(p_val, 3), nsmall = 3))
  # } else if (p_val < 0.05) {
  #  p_text <- paste0("p = ", format(round(p_val, 2), nsmall = 2))
  # } else {
  #  p_text <- paste0("p = ", format(round(p_val, 2), nsmall = 2))
  # }

  # Format correlation coefficient (rounded to 2 decimal places)
  if (is.na(p_val)) {
    p_text = "NA"
    corr_coef <- NA
  } else {
    p_text = format.pval(p_val,digits = 1)
  }
  corr_text <- paste0("r = ", round(corr_coef, 2), ", ", p_text)
  
  # Create plot
  p <- ggplot(data, aes_string(x = x_var, y = y_var, col = col_var)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    annotate("label", x = Inf, y = -Inf, label = corr_text,
             hjust = 1.1, vjust = -0.5) +
    theme_bw() +
    labs(x = x_var, y = y_var)
  p
}

plot_filter_distribution = function(merged_df, colors, log_scale = TRUE) {
  p = merged_df %>% 
    dplyr::group_by(FILTER) %>% 
    dplyr::summarise(n=n()) %>% 
    na.omit() %>% 
    ggplot(mapping = aes(x=reorder(FILTER, +n), y=n, fill = FILTER)) +
    scale_fill_manual(values = colors) + 
    geom_col() +
    theme_bw() +
    coord_flip() +
    labs(x = "FILTER", y="count") + 
    theme(legend.position='None')
  
  if (log_scale) {
    p <- p +
      scale_y_continuous(trans = "log10")
  }
  
  p
}

# Function to calculate Precision-Recall Curve
plot_precision_recall <- function(df) {
  pr <- df %>%
    dplyr::arrange(desc(VAF_caller)) %>%
    dplyr::mutate(
      TP = cumsum(true_positive),
      FP = cumsum(false_positive),
      FN = sum(true_positive) - TP,
      precision = TP / (TP + FP),
      recall = TP / (TP + FN)
    )
  
  # Find the closest point to VAF threshold = 0.1
  vaf_threshold_point <- df %>%
    dplyr::filter(!is.na(VAF_caller)) %>%
    dplyr::slice(which.min(abs(VAF_caller - 0.1)))  # Find closest VAF to 0.1
  
  # Find the corresponding precision and recall at this point
  threshold_recall <- pr$recall[which.min(abs(df$VAF_caller - 0.1))]
  threshold_precision <- pr$precision[which.min(abs(df$VAF_caller - 0.1))]
  
  ggplot(pr, aes(x = recall, y = precision)) +
    geom_line(color = "blue") +
    geom_point(aes(x = threshold_recall, y = threshold_precision), 
               color = "red", size = 3) +  # Add red point at VAF = 0.1
    annotate("text", x = threshold_recall, y = threshold_precision, 
             label = "VAF = 0.1", vjust = -1, color = "red", size = 4) +
    labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") +
    theme_bw()
}

plot_roc_curve <- function(df, vaf_thresh = 0.1) {
  roc_data <- df %>%
    arrange(desc(VAF_caller)) %>%
    mutate(
      TP = cumsum(true_positive),
      FP = cumsum(false_positive),
      TN = sum(true_negative),
      FN = sum(false_negative),
      TPR = TP / (TP + FN),  # Sensitivity / Recall
      FPR = FP / (FP + TN)   # False Positive Rate
    )
  
  # Find the closest point to VAF threshold = vaf_thresh
  vaf_threshold_point <- df %>%
    dplyr::filter(!is.na(VAF_caller)) %>%
    dplyr::slice(which.min(abs(VAF_caller - vaf_thresh)))
  
  # Get corresponding TPR and FPR at this point
  threshold_tpr <- roc_data$TPR[which.min(abs(df$VAF_caller - vaf_thresh))]
  threshold_fpr <- roc_data$FPR[which.min(abs(df$VAF_caller - vaf_thresh))]
  
  ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Diagonal line for random classifier
    geom_point(aes(x = threshold_fpr, y = threshold_tpr), 
               color = "red", size = 3) +  # Highlight VAF = 0.1
    annotate("text", x = threshold_fpr, y = threshold_tpr, 
             label = paste0("VAF = ", vaf_thresh), vjust = -1, color = "red", size = 4) +
    labs(title = "ROC Curve", x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)") +
    theme_bw()
}

# Function to plot Histogram of VAF Differences

plot_vaf_difference <- function(df) {
  df <- df %>%
    mutate(VAF_diff = VAF_caller - VAF_truth)
  
  ggplot(df, mapping = aes(x=FILTER, y=VAF_diff)) +
    geom_boxplot() +
    labs(title = "Boxplot of VAF differences", y = "VAF Difference (Caller - Truth)", x = "") +
    theme_bw()
}

# Function to plot Histogram of depth Differences
plot_cov_difference <- function(df) {
  df <- df %>%
    dplyr::mutate(DP_diff = DP_caller - DP_truth)
  
  ggplot(df, mapping = aes(x=FILTER, y=DP_diff)) +
    geom_boxplot() +
    labs(title = "Boxplot of DP differences", y = "DP Difference (Caller - Truth)", x = "") +
    theme_bw()
}

# Function to plot Read Depth vs. SNP Calling Performance
plot_depth_vs_performance <- function(df) {
  df %>%
    group_by(DP_caller) %>%
    summarize(precision = mean(true_positive / (true_positive + false_positive), na.rm = TRUE),
              recall = mean(true_positive / (true_positive + false_negative), na.rm = TRUE)) %>%
    pivot_longer(cols = c(precision, recall), names_to = "metric", values_to = "value") %>%
    ggplot(aes(x = DP_caller, y = value, color = metric)) +
    geom_line() +
    labs(title = "Read Depth vs. Precision/Recall", x = "Read Depth", y = "Performance Metric") +
    theme_bw()
}

# Function to plot Chromosomal Distribution of Errors
plot_chromosomal_errors <- function(df) {
  df %>%
    group_by(chr_caller) %>%
    summarize(false_positive = sum(false_positive), false_negative = sum(false_negative)) %>%
    pivot_longer(cols = c(false_positive, false_negative), names_to = "error_type", values_to = "count") %>%
    ggplot(aes(x = chr_caller, y = count, fill = error_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Chromosomal Distribution of Errors", x = "Chromosome", y = "Count") +
    theme_bw()
}

# Function to plot Depth Comparison
plot_depth_comparison <- function(df) {
  df <- df %>%
    filter(!is.na(DP_caller) & !is.na(DP_truth))  # Remove NA depths
  
  ggplot(df, aes(x = DP_truth, y = DP_caller)) +
    geom_point(alpha = 0.5, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Comparison of Read Depth (Caller vs. Ground Truth)", 
         x = "Ground Truth Depth", y = "Caller Depth") +
    theme_bw()
}

# Function to plot races coverage histogram
plot_races_coverage <- function(seq_res_long) {
  seq_res_long %>% 
    ggplot(mapping = aes(x=DP)) +
    geom_histogram(binwidth = 1) +
    theme_bw() +
    geom_vline(aes(xintercept = median(DP)), color = "indianred") +
    labs(title = "rRACES coverage",
         subtitle = paste0(nrow(seq_res_long), " mutations; Median coverage = ",
                           median(seq_res_long$DP)),
         color = "")
}

# Function to plot races VAF histogram
plot_races_vaf <- function(seq_res_long) {
  seq_res_long %>% 
    ggplot(mapping = aes(x=VAF)) +
    xlim(c(0,1)) +
    geom_histogram(binwidth = 0.01) +
    theme_bw() +
    labs(title = "rRACES VAF",
         subtitle = paste0(nrow(seq_res_long), " mutations"),
         color = "")
}

# Function to plot caller coverage histogram
plot_caller_coverage <- function(caller_res_pass, sample_info) {
  caller_res_pass %>% 
    ggplot(mapping = aes(x=DP)) +
    geom_histogram(binwidth = 1) +
    theme_bw() +
    geom_vline(aes(xintercept = median(DP)), color = "indianred") +
    labs(title = paste0(sample_info$caller_name, " coverage"),
         subtitle = paste0(nrow(caller_res_pass), 
                           " PASS mutations; Median coverage = ", 
                           median(caller_res_pass$DP)),
         color = "")
}

# Function to plot caller VAF histogram
plot_caller_vaf <- function(caller_res_pass, sample_info) {
  caller_res_pass %>% 
    ggplot(mapping = aes(x=VAF)) +
    geom_histogram(binwidth = 0.01) +
    theme_bw() +
    xlim(c(0,1)) +
    labs(title = paste0(sample_info$caller_name, " VAF"),
         subtitle = paste0(nrow(caller_res_pass), " PASS mutations"),
         color = "")
}

# Function to plot false negative VAF distribution
plot_false_negative_vaf_dist <- function(merged_df0) {
  merged_df0 = merged_df0 %>% 
    dplyr::filter(false_negative)
  median_VAF0 = median(merged_df0$VAF_truth)
  
  merged_df0 %>% 
    dplyr::filter(false_negative) %>% 
    ggplot(mapping = aes(x=VAF_truth)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept=median_VAF0, col="indianred") +
    scale_x_continuous(transform = "log10") +
    theme_bw() +
    labs(x="VAF races", y="Count", title="False negatives' Variant allele frequency")
}

# Function to plot VAF scatter with correlation for all mutations
plot_vaf_scatter_all <- function(merged_df, colors, sample_info) {
  merged_df = merged_df %>% dplyr::filter(positive_truth | positive_call)
  p <- plot_scatter_with_corr(merged_df, "VAF_truth", "VAF_caller") +
    ggtitle("VAF correlation", subtitle = paste0(nrow(merged_df), " total mutations")) +
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    theme(legend.position = "bottom") +
    labs(
      x = "VAF races",
      y = paste0("VAF ", sample_info$caller_name)
    ) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  p_with_marginal <- ggMarginal(p, type = "boxplot", groupFill = TRUE, groupColour = TRUE)
  ggplotify::as.ggplot(p_with_marginal)
}

# Function to plot DP scatter with correlation for all mutations
plot_dp_scatter_all <- function(merged_df, colors, sample_info) {
  merged_df = merged_df %>% dplyr::filter(positive_truth | positive_call)
  p <- plot_scatter_with_corr(merged_df, "DP_truth", "DP_caller") +
    ggtitle("DP correlation", subtitle = paste0(nrow(merged_df), " total mutations")) +
    theme(legend.position = "bottom") +
    labs(
      x = "DP races",
      y = paste0("DP ", sample_info$caller_name)
    ) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  p_with_marginal <- ggMarginal(
    p,
    type = "boxplot",
    groupFill = TRUE, 
    groupColour = TRUE
  )
  
  ggplotify::as.ggplot(p_with_marginal)
}

# Function to plot VAF scatter with correlation for PASS mutations
plot_vaf_scatter_pass <- function(merged_df, colors, sample_info) {
  merged_df = merged_df %>% dplyr::filter(positive_truth | positive_call)
  p <- plot_scatter_with_corr(merged_df, "VAF_truth", "VAF_caller") +
    ggtitle("VAF correlation", subtitle = paste0(nrow(merged_df), " total mutations")) +
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    theme(legend.position = "bottom") +
    labs(
      x = "VAF races",
      y = paste0("VAF ", sample_info$caller_name)
    ) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  p_with_marginal <- ggMarginal(p, type = "boxplot", groupFill = TRUE, groupColour = TRUE)
  ggplotify::as.ggplot(p_with_marginal)
}

# Function to plot DP scatter with correlation for PASS mutations
plot_dp_scatter_pass <- function(merged_df, colors, sample_info) {
  merged_df = merged_df %>% dplyr::filter(positive_truth | positive_call)
  p <- plot_scatter_with_corr(merged_df, "DP_truth", "DP_caller") +
    ggtitle("DP correlation", subtitle = paste0(nrow(merged_df), " total mutations")) +
    theme(legend.position = "bottom") +
    labs(
      x = "DP races",
      y = paste0("DP ", sample_info$caller_name)
    ) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  p_with_marginal <- ggMarginal(
    p,
    type = "boxplot",
    groupFill = TRUE, 
    groupColour = TRUE
  )
  
  ggplotify::as.ggplot(p_with_marginal)
}

# Function to plot metrics comparison
plot_metrics_comparison <- function(metrics) {
  metrics %>% 
    ggplot(mapping = aes(x=name, y=value, fill=Mutations)) +
    geom_col(position = "dodge") +
    theme_bw() +
    ylim(c(0,1)) +
    labs(title = "Performance Metrics Comparison", 
         x = "Metric", 
         y = "Value")
}

# Helper function to compute metrics from vectors
compute_metrics_from_vectors <- function(y_true, y_pred) {
  # Create confusion matrix
  cm <- table(Actual = y_true, Predicted = y_pred)
  
  # Calculate metrics
  TP <- sum(y_true == 1 & y_pred == 1)
  FP <- sum(y_true == 0 & y_pred == 1)
  TN <- sum(y_true == 0 & y_pred == 0)
  FN <- sum(y_true == 1 & y_pred == 0)
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  sensitivity <- TP / (TP + FN)
  precision <- TP / (TP + FP)
  recall <- sensitivity
  f1_score <- 2 * (precision * recall) / (precision + recall)
  false_positive_rate = FP / (FP + TN)
  
  metrics <- dplyr::tibble(
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score,
    FPR = false_positive_rate
  )
  
  return(metrics)
}

# Plot metric over VAF threshold
plot_metric_over_VAF_threshold = function(seq_res_long, caller_res, only_pass, VAF_spectrum=seq(0,.1, by=.005)) {
  if (only_pass) {
    caller_res = caller_res %>% dplyr::filter(FILTER == "PASS")
  }

  # dataframe with metrics (dfm)  
  dfm = lapply(VAF_spectrum, function(min_vaf) {
    merged_df = merge_datasets(caller_res, seq_res_long, min_vaf)
    y_true <- as.numeric(factor((merged_df$VAF_truth > min_vaf), levels=c(FALSE, TRUE))) - 1
    y_pred <- as.numeric(factor((merged_df$VAF_caller > min_vaf), levels=c(FALSE, TRUE))) - 1
    metrics = compute_metrics_from_vectors(y_true, y_pred) 
    dplyr::tibble(min_vaf=min_vaf, value = as.numeric(metrics), metric = colnames(metrics))  
  }) %>% do.call("bind_rows", .)
  
  metric_colors = c(
    "#1f77b4",    # Blue
    "#ff7f0e",       # Orange
    "#2ca02c",     # Green
    "#d62728",  # Red
    "#9467bd"      # Purple
  )
  
  dfm %>% 
    ggplot(mapping = aes(x=min_vaf, y=value, col=metric)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    labs(x = "VAF threshold", y="Value") +
    scale_color_manual(values = metric_colors) +
    #theme(legend.position = "none") +
    ylim(c(0,1))
}

plot_flow_of_calls = function(merged_df) {
  merged_df$FILTER[is.na(merged_df$FILTER)] = "Not called"
  df = merged_df %>% 
    dplyr::select(FILTER, positive_truth, positive_call) %>% 
    dplyr::filter(positive_truth | positive_call) %>% 
    dplyr::mutate(
      Pre_filtering = dplyr::case_when(
        positive_truth & positive_call ~ "shared",
        positive_truth ~ "races private",
        positive_call ~ "caller private",
        .default = "N/A"
      ),
      Post_filtering = dplyr::case_when(
        (FILTER != "PASS") & (!positive_truth) ~ "Correclty filtered out",
        (FILTER != "PASS") & (positive_truth) ~ "races private",
        positive_truth & positive_call ~ "shared",
        positive_truth ~ "races private",
        positive_call ~ "caller private",
        .default = "N/A"
      )
    ) %>% 
    dplyr::select(FILTER, Pre_filtering, Post_filtering) %>% 
    dplyr::group_by(FILTER, Pre_filtering, Post_filtering) %>% 
    dplyr::summarise(Freq = n())
  
  df %>% 
    #dplyr::mutate(FILTER = ifelse(Pre_filtering %in% c("races private", "shared"), "Correct", "Wrong")) %>% 
    ggplot(mapping = aes(y=Freq, axis1=Pre_filtering, axis2=Post_filtering)) +
    geom_flow(aes(fill=FILTER), curve_type = "quintic") +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Pre_filtering", "Post_filtering"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    theme_bw()  
  
  
  
  
  df = merged_df %>% 
    dplyr::select(FILTER, positive_truth, positive_call) %>% 
    dplyr::filter(positive_truth | positive_call) %>% 
    dplyr::mutate(
      Pre_filtering = dplyr::case_when(
        positive_truth & positive_call ~ "TP",
        positive_truth ~ "FN",
        positive_call ~ "FP",
        .default = "TN"
      ),
      Post_filtering = dplyr::case_when(
        (FILTER != "PASS") & (!positive_truth) ~ "TN",
        (FILTER != "PASS") & (positive_truth) ~ "FN",
        positive_truth & positive_call ~ "TP",
        positive_truth ~ "FN",
        positive_call ~ "FP",
        .default = "N/A"
      )
    ) %>% 
    dplyr::select(FILTER, Pre_filtering, Post_filtering) %>% 
    dplyr::group_by(FILTER, Pre_filtering, Post_filtering) %>% 
    dplyr::summarise(Freq = n())
  
  df %>% 
    #dplyr::mutate(FILTER = ifelse(Pre_filtering %in% c("races private", "shared"), "Correct", "Wrong")) %>% 
    ggplot(mapping = aes(y=Freq, axis1=Pre_filtering, axis2=Post_filtering)) +
    geom_flow(aes(fill=FILTER), curve_type = "quintic") +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Pre_filtering", "Post_filtering"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    theme_bw()  
}


plot_ecdf_comparison <- function(vec1, vec2, label1, label2, x_label) {
  library(ggplot2)
  library(dplyr)
  
  # Combine vectors into a data frame
  df <- bind_rows(
    data.frame(value = vec1, group = label1),
    data.frame(value = vec2, group = label2)
  )
  
  # Perform the KS test
  ks_result <- ks.test(vec1, vec2)
  
  # Plot ECDF
  p <- ggplot(df, aes(x = value, color = group)) +
    stat_ecdf() +
    labs(
      title = "ECDF Comparison",
      subtitle = sprintf("KS test p-value: %.4g", ks_result$p.value),
      x = x_label,
      y = "ECDF",
      color = ""
    ) +
    theme_minimal()
  
  return(p)
}

plot_density_comparison <- function(vec1, vec2, label1, label2, x_label) {
  library(ggplot2)
  library(dplyr)
  
  # Combine vectors into a data frame
  df <- bind_rows(
    data.frame(value = vec1, group = label1),
    data.frame(value = vec2, group = label2)
  )
  
  # Perform the KS test
  ks_result <- ks.test(vec1, vec2)
  
  # Sample sizes
  # Format sample sizes with apostrophes
  n1 <- format(length(vec1), big.mark = "'")
  n2 <- format(length(vec2), big.mark = "'")
  
  # Build subtitle with formatted sample sizes
  subtitle_text <- sprintf(
    "%s (n = %s), %s (n = %s)\nKS test p-value: %.4g",
    label1, n1, label2, n2, ks_result$p.value
  )
  
  # Plot densities
  p <- ggplot(df, aes(x = value, fill = group)) +
    geom_density(alpha = 0.5) +
    labs(
      title = paste0(x_label, " comparison"),
      subtitle = subtitle_text,
      x = x_label,
      y = "Density",
      fill = ""
    ) +
    theme_bw()
  
  p
}



plot_density_comparison_multi <- function(vec_list, labels, colors, x_label) {
  library(ggplot2)
  library(dplyr)
  library(purrr)
  
  if (length(vec_list) != length(labels) || length(labels) != length(colors)) {
    stop("vec_list, labels, and colors must have the same length.")
  }
  
  # Combine into one data frame
  df <- map2_dfr(vec_list, labels, ~ data.frame(value = .x, group = .y))
  
  # Pairwise KS test (just first vs. rest)
  ks_tests <- map2(vec_list[-1], labels[-1], ~ {
    ks.test(vec_list[[1]], .x)
  })
  p_values <- map_chr(ks_tests, ~ sprintf("%.4g", .x$p.value))
  comp_labels <- paste0(labels[1], " vs ", labels[-1], ": p = ", p_values)
  
  # Build subtitle with n values and p-values
  ns <- map_chr(vec_list, ~ format(length(.x), big.mark = "'"))
  group_info <- paste0(labels, " (n = ", ns, ")")
  subtitle_text <- paste(
    paste(group_info, collapse = ", "),
    paste(comp_labels, collapse = " | "),
    sep = "\n"
  )
  
  subtitle_text <- paste(
    paste(group_info, collapse = ", ")
  )
  
  # Plot
  p <- ggplot(df, aes(x = value, fill = group)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = colors) +
    labs(
      title = paste0(x_label, " Comparison"),
      subtitle = subtitle_text,
      x = x_label,
      y = "Density",
      fill = ""
    ) +
    theme_bw()
  
  p
}


plot_ecdf_comparison_multi <- function(vec_list, labels, colors, x_label) {
  library(ggplot2)
  library(dplyr)
  library(purrr)
  
  if (length(vec_list) != length(labels) || length(labels) != length(colors)) {
    stop("vec_list, labels, and colors must have the same length.")
  }
  
  # Combine into one data frame
  df <- map2_dfr(vec_list, labels, ~ data.frame(value = .x, group = .y))
  
  # Pairwise KS test (first vs rest)
  ks_tests <- map2(vec_list[-1], labels[-1], ~ {
    ks.test(vec_list[[1]], .x)
  })
  p_values <- map_chr(ks_tests, ~ sprintf("%.4g", .x$p.value))
  comp_labels <- paste0(labels[1], " vs ", labels[-1], ": p = ", p_values)
  
  # Subtitle with sample sizes and p-values
  ns <- map_chr(vec_list, ~ format(length(.x), big.mark = "'"))
  group_info <- paste0(labels, " (n = ", ns, ")")
  subtitle_text <- paste(
    paste(group_info, collapse = ", "),
    paste(comp_labels, collapse = " | "),
    sep = "\n"
  )
  
  subtitle_text <- paste(
    paste(group_info, collapse = ", ")
  )
  
  # Plot
  p <- ggplot(df, aes(x = value, color = group)) +
    stat_ecdf(geom = "step", size = 1) +
    scale_color_manual(values = colors) +
    labs(
      title = paste0(x_label, " ECDF Comparison"),
      subtitle = subtitle_text,
      x = x_label,
      y = "ECDF",
      color = ""
    ) +
    theme_bw()
  
  p
}


plot_precision_recall_vaf <- function(vaf_analysis_results, 
                                      title = "Precision and Recall Across VAF Bins",
                                      text_size = 3.5,
                                      point_size = 1.5,
                                      line_size = 1) {
  
  # Extract performance table from analysis results
  performance_data <- vaf_analysis_results$performance_table
  
  # Check if required columns exist
  required_cols <- c("VAF_bin", "precision", "sensitivity")
  if (!all(required_cols %in% colnames(performance_data))) {
    stop("Required columns missing from performance table: ", 
         paste(setdiff(required_cols, colnames(performance_data)), collapse = ", "))
  }
  
  # Prepare data for plotting (rename sensitivity to recall for clarity)
  plot_data <- performance_data %>%
    dplyr::select(VAF_bin, precision, recall = sensitivity) %>%
    dplyr::filter(!is.na(VAF_bin)) %>%
    tidyr::pivot_longer(cols = c(precision, recall), 
                        names_to = "metric", 
                        values_to = "value") %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("precision", "recall"), 
                      labels = c("Precision", "Recall")),
      value_pct = round(value * 100, 1),
      label_text = paste0(value_pct, "%")
    )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = VAF_bin, y = value, color = metric, group = metric)) +
    # Add lines connecting points
    geom_line(size = line_size, alpha = 0.8) +
    # Add points
    geom_point(size = point_size, alpha = 0.9) +
    # Add text labels showing percentages
    geom_text(aes(label = label_text), 
              vjust = -0.5, 
              size = text_size, 
              show.legend = FALSE,
              fontface = "bold") +
    # Set y-axis to 0-1 range
    scale_y_continuous(limits = c(0, 1.05), 
                       breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(accuracy = 1)) +
    # Color scheme
    scale_color_manual(values = c("Precision" = "#E31A1C", "Recall" = "#1F78B4"),
                       name = "Metric") +
    # Labels and title
    labs(
      title = title,
      x = "VAF Bin",
      y = "Performance",
      caption = paste0("VAF tolerance: ", vaf_analysis_results$vaf_tolerance_pct, "%; ",
                       "Min VAF threshold: ", vaf_analysis_results$min_vaf_threshold)
    ) +
    # Theme
    theme_bw()
  p
  
  # Add overall metrics as annotation if available
  if (!is.null(vaf_analysis_results$overall_metrics)) {
    overall <- vaf_analysis_results$overall_metrics
    if (!is.na(overall$precision) && !is.na(overall$sensitivity)) {
      annotation_text <- paste0(
        "Overall: Precision = ", round(overall$precision * 100, 1), "%, ",
        "Recall = ", round(overall$sensitivity * 100, 1), "%"
      )
      
      p <- p + labs(subtitle = annotation_text)
    }
  }
  
  return(p)
}

# Create combined sensitivity and FPR plot
plot_across_vaf <- function(results) {
  plot_data <- results$performance_table %>%
    select(VAF_bin, sensitivity_pct, precision_pct) %>%
    #filter(!is.na(sensitivity_pct) | !is.na(fpr_pct)) %>%  # Remove bins with no data
    tidyr::pivot_longer(cols = c(sensitivity_pct, precision_pct), 
                        names_to = "metric", values_to = "percentage")
  
  ggplot(plot_data, aes(x = VAF_bin, y = percentage, fill = metric)) +
    geom_col(position = "dodge", alpha = 0.7) +
    geom_text(aes(label = ifelse(!is.na(percentage), paste0(round(percentage, 1), "%"), "NA")), 
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    labs(
      title = "Variant Caller Performance by VAF Range",
      x = "VAF Range (Truth)",
      y = "Percentage (%)",
      fill = "Metric"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("sensitivity_pct" = "steelblue", "precision_pct" = "coral"),
                      labels = c("Recall", "Precision")) +
    ylim(0, max(plot_data$percentage, na.rm = TRUE) * 1.1) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100))
}

# Create separate FPR plot with improved TN definition
plot_fpr_by_vaf <- function(results) {
  plot_data <- results$performance_table %>%
    filter(!is.na(fpr_pct))
  
  ggplot(plot_data, aes(x = VAF_bin, y = fpr_pct)) +
    geom_col(fill = "coral", alpha = 0.7) +
    geom_text(aes(label = paste0(round(fpr_pct, 2), "%")), vjust = -0.5) +
    labs(
      title = "False Positive Rate by VAF Range",
      x = "VAF Range",
      y = "False Positive Rate (%)",
      subtitle = "True Negatives = Ground truth variants below minimum VAF threshold"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, max(plot_data$fpr_pct, na.rm = TRUE) * 1.1)
}

# Create specificity plot (complement of FPR)
plot_specificity_by_vaf <- function(results) {
  plot_data <- results$performance_table %>%
    filter(!is.na(specificity_pct))
  
  ggplot(plot_data, aes(x = VAF_bin, y = specificity_pct)) +
    geom_col(fill = "lightgreen", alpha = 0.7) +
    geom_text(aes(label = paste0(round(specificity_pct, 1), "%")), vjust = -0.5) +
    labs(
      title = "Specificity by VAF Range",
      x = "VAF Range",
      y = "Specificity (%)",
      subtitle = "Ability to correctly identify true negatives (low VAF variants)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 100)
}

# Create VAF accuracy plot
plot_vaf_accuracy <- function(results) {
  accuracy_data <- results$performance_table %>%
    select(VAF_bin, mean_abs_error, median_abs_error, n_variants) %>%
    filter(!is.na(mean_abs_error) & n_variants > 0) %>%
    tidyr::pivot_longer(cols = c(mean_abs_error, median_abs_error), 
                        names_to = "metric", values_to = "error")
  
  ggplot(accuracy_data, aes(x = VAF_bin, y = error, fill = metric)) +
    geom_col(position = "dodge", alpha = 0.7) +
    labs(
      title = "VAF Calling Accuracy by VAF Range (True Positives Only)",
      x = "VAF Range (Truth)", 
      y = "Absolute Error",
      fill = "Error Metric"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("mean_abs_error" = "coral", "median_abs_error" = "lightblue"),
                      labels = c("Mean Abs Error", "Median Abs Error"))
}



plot_multicaller_precision_recall_vaf <- function(performance_data, 
                                                  title = "Precision and Recall Across VAF Bins",
                                                  text_size = 3.5,
                                                  point_size = 1.5,
                                                  line_size = 1) {
  
  # Check if required columns exist
  required_cols <- c("VAF_bin", "precision", "sensitivity", "caller")
  if (!all(required_cols %in% colnames(performance_data))) {
    stop("Required columns missing from performance table: ", 
         paste(setdiff(required_cols, colnames(performance_data)), collapse = ", "))
  }
  
  # Prepare data for plotting (rename sensitivity to recall for clarity)
  plot_data <- performance_data %>%
    dplyr::select(VAF_bin, precision, recall = sensitivity, caller) %>%
    dplyr::filter(!is.na(VAF_bin)) %>%
    tidyr::pivot_longer(cols = c(precision, recall), 
                        names_to = "metric", 
                        values_to = "value") %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("precision", "recall"), 
                      labels = c("Precision", "Recall")),
      value_pct = round(value * 100, 1),
      label_text = paste0(value_pct, "%")
    )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = VAF_bin, y = value, color = caller, linetype = metric, group = interaction(caller, metric))) +
    # Add lines connecting points
    geom_line(size = line_size, alpha = 0.8) +
    # Add points
    geom_point(size = point_size, alpha = 0.9) +
    # Add text labels showing percentages
    geom_text(aes(label = label_text), 
              vjust = -0.5, 
              size = text_size, 
              show.legend = FALSE,
              fontface = "bold") +
    # Set y-axis to 0-1 range
    scale_y_continuous(limits = c(0, 1.05), 
                       breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(accuracy = 1)) +
    # Linetype for metrics
    scale_linetype_manual(values = c("Precision" = "solid", "Recall" = "dashed"),
                          name = "Metric") +
    # Color scheme for callers (you may want to adjust these colors based on your callers)
    scale_color_manual(name = "Caller", values = method_colors) +
    # Labels and title
    labs(
      title = title,
      x = "VAF Bin",
      y = "Performance"
    ) +
    # Theme
    theme_bw()
  p
  
  return(p)
}


plot_multicaller_rmse_vaf <- function(performance_data, 
                                      title = "RMSE Across VAF Bins",
                                      text_size = 3.5,
                                      point_size = 1.5,
                                      line_size = 1) {
  
  # Check if required columns exist
  required_cols <- c("VAF_bin", "rmse", "caller")
  if (!all(required_cols %in% colnames(performance_data))) {
    stop("Required columns missing from performance table: ", 
         paste(setdiff(required_cols, colnames(performance_data)), collapse = ", "))
  }
  
  # Prepare data for plotting
  plot_data <- performance_data %>%
    dplyr::select(VAF_bin, rmse, caller) %>%
    dplyr::filter(!is.na(VAF_bin), !is.na(rmse)) %>%
    dplyr::mutate(
      rmse_rounded = round(rmse, 3),
      label_text = as.character(rmse_rounded)
    )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = VAF_bin, y = rmse, color = caller, group = caller)) +
    # Add lines connecting points
    geom_line(size = line_size, alpha = 0.8) +
    # Add points
    geom_point(size = point_size, alpha = 0.9) +
    # Add text labels showing RMSE values
    # geom_text(aes(label = label_text), 
    #           vjust = -0.5, 
    #           size = text_size, 
    #           show.legend = FALSE,
    #           fontface = "bold") +
    # Set y-axis to start from 0
    scale_y_continuous(limits = c(0, max(plot_data$rmse, na.rm = TRUE) * 1.1)) +
    # Color scheme for callers
    scale_color_manual(name = "Caller", values = method_colors) +
    # Labels and title
    labs(
      title = title,
      x = "VAF Bin",
      y = "RMSE"
    ) +
    # Theme
    theme_bw()
  
  p
}

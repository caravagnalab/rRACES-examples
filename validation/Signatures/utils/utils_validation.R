evaluate_signatures <- function(true_signatures, predicted_signatures) {
  # Union of all signatures involved
  all_signatures <- union(true_signatures, predicted_signatures)
  
  # Binary vectors
  true_labels <- ifelse(all_signatures %in% true_signatures, 1, 0)
  pred_labels <- ifelse(all_signatures %in% predicted_signatures, 1, 0)
  
  # Confusion matrix
  cm <- table(True = true_labels, Predicted = pred_labels)
  print("Confusion Matrix:")
  print(cm)
  
  # Compute metrics
  TP <- sum(true_labels == 1 & pred_labels == 1)
  FP <- sum(true_labels == 0 & pred_labels == 1)
  FN <- sum(true_labels == 1 & pred_labels == 0)
  TN <- sum(true_labels == 0 & pred_labels == 0)
  
  accuracy <- (TP + TN) / (TP + FP + FN + TN)
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  f1_score <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA
  }
  
  return(list(
    ConfusionMatrix = cm,
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN
  ))
}


plot_confusion_matrix <- function(conf_matrix, tool_name) {
  conf_df <- as.data.frame(as.table(conf_matrix))
  colnames(conf_df) <- c("True", "Predicted", "Count")

  conf_df$True <- factor(conf_df$True, levels = c("1", "0"))
  conf_df$Predicted <- factor(conf_df$Predicted, levels = c("1", "0"))

  conf_df$Label <- with(conf_df, ifelse(Predicted == "1" & True == "1", "TP",
                                        ifelse(Predicted == "0" & True == "0", "TN",
                                               ifelse(Predicted == "1" & True == "0", "FP", "FN"))))

  conf_df$Group <- ifelse(conf_df$Label %in% c("TP", "TN"), "Correct", "Incorrect")
  color_map <- c("Correct" = "#91D1C2B2", "Incorrect" = "#F39B7FB2")

  ggplot(conf_df, aes(x = Predicted, y = True, fill = Group)) +
    geom_tile(color = "black") +
    geom_text(aes(label = paste0(Label, "\n", Count)), size = 6) +
    scale_fill_manual(values = color_map, name = "Classification") +
    scale_x_discrete(labels = c("0" = "Negative", "1" = "Positive")) +
    scale_y_discrete(labels = c("1" = "Positive", "0" = "Negative")) +
    labs(title = paste(tool_name),
         x = "Predicted Label",
         y = "True Label") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "right",
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 11))
}


align_columns <- function(df, target_columns) {
  df <- df[, intersect(colnames(df), target_columns), drop = FALSE]
  # Add missing columns as zero
  missing_cols <- setdiff(target_columns, colnames(df))
  for (col in missing_cols) {
    df[[col]] <- 0
  }

  df <- df[, target_columns, drop = FALSE]
  return(df)
}


cosine_similarity_exposures <- function(vec1, vec2) {
  if (all(vec1 == 0) || all(vec2 == 0)) return(NA)
  return(sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2))))
}


compute_sample_cosine_similarity <- function(mat1, mat2) {
  sapply(seq_len(nrow(mat1)), function(i) {
    cosine_similarity_exposures(mat1[i, ], mat2[i, ])
  })
}


compute_similarity_matrix <- function(inferred_mat, ground_truth_mat) {
  similarity_matrix <- matrix(NA,
                              nrow = ncol(inferred_mat),
                              ncol = ncol(ground_truth_mat),
                              dimnames = list(colnames(inferred_mat), colnames(ground_truth_mat)))

  for (i in seq_len(ncol(inferred_mat))) {
    for (j in seq_len(ncol(ground_truth_mat))) {
      similarity_matrix[i, j] <- cosine_similarity_exposures(inferred_mat[, i], ground_truth_mat[, j])
    }
  }

  return(similarity_matrix)
}


compute_mse <- function(inferred_mat, ground_truth_mat) {
  inferred_mat <- as.matrix(inferred_mat)
  ground_truth_mat <- as.matrix(ground_truth_mat)

  # Align columns
  common_cols <- intersect(colnames(inferred_mat), colnames(ground_truth_mat))
  inferred_mat <- inferred_mat[, common_cols, drop = FALSE]
  ground_truth_mat <- ground_truth_mat[, common_cols, drop = FALSE]

  # Compute squared error matrix
  mse_matrix <- (inferred_mat - ground_truth_mat)^2

  # MSE per signature (column-wise mean)
  mse_per_signature <- colMeans(mse_matrix, na.rm = TRUE)
  names(mse_per_signature) <- common_cols  # Ensure names

  mse_overall <- mean(mse_matrix, na.rm = TRUE)

  return(list(
    per_signature = mse_per_signature,
    overall = mse_overall
  ))
}


validate_exposures <- function(process_exposures, sparsesig_out, sigprofiler_out) {
  sparsesig_exp <- align_columns(sparsesig_out, colnames(process_exposures))
  sigprofiler_exp <- align_columns(sigprofiler_out, colnames(process_exposures))

  sparsesig_exp <- as.matrix(sparsesig_exp)
  sigprofiler_exp <- as.matrix(sigprofiler_exp)
  process_exposures <- as.matrix(process_exposures)

  # Cosine similarity per sample
  similarity_sparsesig <- compute_sample_cosine_similarity(sparsesig_exp, process_exposures)
  similarity_sigprofiler <- compute_sample_cosine_similarity(sigprofiler_exp, process_exposures)

  # MSE overall per sample
  mse_sparsesig <- mean(rowMeans((sparsesig_exp - process_exposures)^2))
  mse_sigprofiler <- mean(rowMeans((sigprofiler_exp - process_exposures)^2))

  # Compute per-signature MSE (across samples)
  mse_per_signature <- function(inferred_mat, true_mat) {
    sapply(seq_len(ncol(inferred_mat)), function(j) {
      mean((inferred_mat[, j] - true_mat[, j])^2)
    })
  }

  mse_sparsesig_per_sig <- mse_per_signature(sparsesig_exp, process_exposures)
  mse_sigprofiler_per_sig <- mse_per_signature(sigprofiler_exp, process_exposures)

  return(list(
    cosine_sparsesig = similarity_sparsesig,
    cosine_sigprofiler = similarity_sigprofiler,
    mse_sparsesig_overall = mse_sparsesig$overall,
    mse_sigprofiler_overall = mse_sigprofiler$overall,
    mse_sparsesig_per_signature = mse_sparsesig$per_signature,
    mse_sigprofiler_per_signature = mse_sigprofiler$per_signature
  ))
}







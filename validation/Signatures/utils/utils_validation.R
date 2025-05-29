generate_confusion_matrix <- function(process_sign, sparsesig_sign, sigprofiler_sign) {
  ground_truth <- colnames(process_sign)[-1]
  sparsesig <- colnames(sparsesig_sign)
  sigprofiler <- colnames(sigprofiler_sign)[-1]
  all_signatures <- unique(c(ground_truth, sparsesig, sigprofiler))
  
  binary_matrix <- function(target, ref) as.numeric(target %in% ref)
  conf_matrix_sparsesig <- table(Predicted = binary_matrix(all_signatures, sparsesig), True = binary_matrix(all_signatures, ground_truth))
  conf_matrix_sigprofiler <- table(Predicted = binary_matrix(all_signatures, sigprofiler), True = binary_matrix(all_signatures, ground_truth))
  
  return(list(conf_matrix_sparsesig = conf_matrix_sparsesig, conf_matrix_sigprofiler = conf_matrix_sigprofiler))
}


plot_confusion_matrix <- function(conf_matrix, tool_name) {
  conf_df <- as.data.frame(as.table(conf_matrix))
  colnames(conf_df) <- c("Predicted", "True", "Count")

  # Assign classification type (TP, TN, FP, FN)
  conf_df$Label <- with(conf_df, ifelse(Predicted == 1 & True == 1, "TP",
                                        ifelse(Predicted == 0 & True == 0, "TN",
                                               ifelse(Predicted == 1 & True == 0, "FP", "FN"))))

  conf_df$Group <- ifelse(conf_df$Label %in% c("TP", "TN"), "True", "False")
  color_map <- c("True" = "#F39B7FB2", "False" = "#91D1C2B2")

  ggplot(conf_df, aes(x = as.factor(Predicted), y = as.factor(True), fill = Group)) +
    geom_tile(color = "black") +
    geom_text(aes(label = paste(Label, "\n", Count)), size = 5) +
    scale_fill_manual(values = color_map, name = "Classification") +
    labs(title = paste(tool_name),
         x = "Predicted",
         y = "True") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 14),
          legend.position = "right")
}


compute_performance_metrics <- function(conf_matrix) {
  TP <- conf_matrix["1", "1"]
  TN <- conf_matrix["0", "0"]
  FP <- conf_matrix["1", "0"]
  FN <- conf_matrix["0", "1"]

  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  Precision <- TP / (TP + FP)
  Recall <- TP / (TP + FN)
  Specificity <- TN / (TN + FP)
  F1_score <- 2 * (Precision * Recall) / (Precision + Recall)

  return(list(
    Accuracy = round(Accuracy, 3),
    Precision = round(Precision, 3),
    Recall = round(Recall, 3),
    Specificity = round(Specificity, 3),
    F1_score = round(F1_score, 3)
  ))
}



align_columns <- function(df, target_columns) {
  df <- df[, colnames(df) %in% target_columns, drop = FALSE]
  missing_cols <- setdiff(target_columns, colnames(df))
  for (col in missing_cols) {
    df[[col]] <- 0
  }
  df <- df[, target_columns, drop = FALSE]

  return(df)
}


normalize_exposures <- function(mat) {
  mat[is.na(mat)] <- 0
  return(sweep(mat, 2, colSums(mat, na.rm = TRUE), FUN = "/"))
}


compute_similarity_matrix <- function(inferred_mat, ground_truth_mat) {
  similarity_matrix <- matrix(NA,
                              nrow = ncol(inferred_mat),
                              ncol = ncol(ground_truth_mat),
                              dimnames = list(colnames(inferred_mat), colnames(ground_truth_mat)))

  # Compute cosine similarity for each signature pairing
  for (i in 1:ncol(inferred_mat)) {
    for (j in 1:ncol(ground_truth_mat)) {
      similarity_matrix[i, j] <- cosine_similarity_exposures(inferred_mat[, i], ground_truth_mat[, j])
    }
  }

  return(similarity_matrix)
}


cosine_similarity_exposures <- function(vec1, vec2) {
  if (any(is.na(vec1)) || any(is.na(vec2))) {
    return(NA)
  }

  if (sum(vec1^2, na.rm = TRUE) == 0 || sum(vec2^2, na.rm = TRUE) == 0) {
    return(NA)
  }

  return(sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2))))
}



compute_mse <- function(inferred_mat, ground_truth_mat) {
  inferred_mat <- as.matrix(sapply(inferred_mat, as.numeric))
  ground_truth_mat <- as.matrix(sapply(ground_truth_mat, as.numeric))

  inferred_mat[is.na(inferred_mat)] <- 0
  ground_truth_mat[is.na(ground_truth_mat)] <- 0

  common_cols <- intersect(colnames(inferred_mat), colnames(ground_truth_mat))
  inferred_mat <- inferred_mat[, common_cols, drop = FALSE]
  ground_truth_mat <- ground_truth_mat[, common_cols, drop = FALSE]

  mse_matrix <- (inferred_mat - ground_truth_mat) ^ 2

  mse_per_signature <- colMeans(mse_matrix, na.rm = TRUE)

  overall_mse <- mean(mse_matrix, na.rm = TRUE)

  return(list(per_signature = mse_per_signature, overall = overall_mse))
}


# Compute exposure similarities #
validate_exposures <- function(process_exposures, sparsesig_out, sigprofiler_out) {

  sparsesig_exp <- align_columns(sparsesig_out, colnames(process_exposures))
  sigprofiler_exp <- align_columns(sigprofiler_out, colnames(process_exposures))

  sparsesig_exp_norm <- normalize_exposures(sparsesig_exp)
  sigprofiler_exp_norm <- normalize_exposures(sigprofiler_exp)
  process_exp_norm <- normalize_exposures(process_exposures)

  similarity_sparsesig <- compute_similarity_matrix(sparsesig_exp_norm, process_exp_norm)
  similarity_sigprofiler <- compute_similarity_matrix(sigprofiler_exp_norm, process_exp_norm)

  mean_cosine_sparsesig <- rowMeans(similarity_sparsesig, na.rm = TRUE)
  mean_cosine_sigprofiler <- rowMeans(similarity_sigprofiler, na.rm = TRUE)

  mse_sparsesig <- compute_mse(sparsesig_exp_norm, process_exp_norm)
  mse_sigprofiler <- compute_mse(sigprofiler_exp_norm, process_exp_norm)

  return(list(
    mean_cosine_sparsesig = mean_cosine_sparsesig,
    mean_cosine_sigprofiler = mean_cosine_sigprofiler,
    mse_sparsesig = mse_sparsesig,
    mse_sigprofiler = mse_sigprofiler
  ))
}

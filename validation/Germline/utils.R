

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

library(caret)

# Function to compute confusion matrix and performance metrics
compute_metrics <- function(actual, predicted) {
  cm <- table(Actual = actual, Predicted = predicted)
  confusion_matrix <- confusionMatrix(as.factor(predicted), as.factor(actual), positive = "1")
  
  metrics <- dplyr::tibble(
    Accuracy = confusion_matrix$overall["Accuracy"],
    Sensitivity = confusion_matrix$byClass["Sensitivity"],
    Specificity = confusion_matrix$byClass["Specificity"],
    Precision = confusion_matrix$byClass["Precision"],
    Recall = confusion_matrix$byClass["Recall"],
    F1_Score = confusion_matrix$byClass["F1"]
  )
  
  return(metrics)
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
    theme_minimal() +
    labs(title = "Confusion Matrix", x = "Predicted", y = "Actual")
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
    geom_point(alpha = 0.8) +
    #geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray90", alpha = 0.3, linewidth = .5) +
    annotate("text", x = Inf, y = -Inf, label = corr_text,
             hjust = 1.1, vjust = -0.5) +
    theme_bw() +
    labs(x = x_var, y = y_var)
  
  p
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


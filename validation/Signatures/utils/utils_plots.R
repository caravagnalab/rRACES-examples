generate_sankey <- function(spn_id, cov = NULL, pur = NULL) {
  df_spn <- sankey_df %>% filter(SPN == spn_id)
  
  # Optional filtering if coverage and purity are provided
  if (!is.null(cov)) df_spn <- df_spn %>% filter(Coverage == cov)
  if (!is.null(pur)) df_spn <- df_spn %>% filter(Purity == pur)
  
  gt_sigs <- unique(df_spn$Signature[df_spn$Method == "GroundTruth"])
  ss_sigs <- unique(df_spn$Signature[df_spn$Method == "SparseSignatures"])
  sp_sigs <- unique(df_spn$Signature[df_spn$Method == "SigProfiler"])
  
  all_sigs <- sort(unique(c(gt_sigs, ss_sigs, sp_sigs)))
  
  presence_df <- data.frame(Signature = all_sigs) %>%
    dplyr::mutate(
      GroundTruth = Signature %in% gt_sigs,
      SparseSignatures = Signature %in% ss_sigs,
      SigProfiler = Signature %in% sp_sigs
    )
  
  long_df <- presence_df %>%
    tidyr::pivot_longer(cols = c(GroundTruth, SparseSignatures, SigProfiler),
                 names_to = "Method",
                 values_to = "Present") %>%
    dplyr::filter(Present) %>%
    dplyr::select(-Present) %>%
    dplyr::mutate(Method = factor(Method, levels = c("GroundTruth", "SparseSignatures", "SigProfiler")),
           alluvium = Signature)
  
  # Build dynamic title
  title_text <- paste0(spn_id)
  if (!is.null(cov)) title_text <- paste0(title_text, ", cov=", cov)
  if (!is.null(pur)) title_text <- paste0(title_text, ", pur=", pur)
  
  # Plot
  ggplot(long_df,
         aes(x = Method, stratum = Signature, alluvium = alluvium,
             fill = Signature, label = Signature)) +
    geom_flow(width = 1/4) +
    geom_stratum(width = 1/3, color = "black") +
    geom_text(stat = "stratum", size = 3) +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, face = "plain")
    ) +
    labs(
      title = title_text,
      x = "Method",
      y = "Signatures"
    )
}


plot_cosine_similarity <- function(data,
                                   x = "Tool",
                                   y = "CosineSimilarity",
                                   color_var = "Sample",
                                   shape_var = "Purity",
                                   facet_row = "Coverage",
                                   facet_col = "SPN",
                                   shape_values = c(16, 17, 15, 18, 3, 4),
                                   title = "Sample-wise Cosine similarity") {

  # Check required columns
  required_cols <- c(x, y, color_var, shape_var, facet_row, facet_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if(length(missing_cols) > 0){
    stop("Data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Convert shape_var to factor
  data <- data %>%
    mutate(!!shape_var := factor(.data[[shape_var]]))

  xq <- sym(x)
  yq <- sym(y)
  colorq <- sym(color_var)
  shapeq <- sym(shape_var)
  facet_row_q <- sym(facet_row)
  facet_col_q <- sym(facet_col)

  # Use as_labeller for Coverage to prepend "Coverage: "
  coverage_labeller <- as_labeller(function(x) paste0("Coverage: ", x))

  ggplot(data, aes(x = !!xq, y = !!yq, color = !!colorq, shape = !!shapeq)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.9, size = 5) +
    facet_grid(rows = vars(!!facet_row_q), cols = vars(!!facet_col_q),
               labeller = labeller(
                 !!facet_col := label_value,
                 !!facet_row := coverage_labeller
               )) +
    scale_shape_manual(values = shape_values) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 16, face = "plain"),
      strip.text.x = element_text(size = 16, face = "bold"),
      strip.text.y = element_text(size = 16, face = "plain")
    ) +
    labs(
      title = title,
      x = x,
      y = y,
      color = color_var,
      shape = shape_var
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3)))
}


plot_mse_per_signature <- function(data,
                                   x = "Signature",
                                   y = "MSE",
                                   color_var = "Tool",
                                   shape_var = "Purity",
                                   facet_row = "Coverage",
                                   facet_col = "SPN",
                                   shape_values = c(16, 17, 15, 18, 3, 4),
                                   title = "MSE per Signature") {

  # Check columns
  required_cols <- c(x, y, color_var, shape_var, facet_row, facet_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if(length(missing_cols) > 0){
    stop("Data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  data <- data %>%
    mutate(!!shape_var := factor(.data[[shape_var]]))

  xq <- sym(x)
  yq <- sym(y)
  colorq <- sym(color_var)
  shapeq <- sym(shape_var)
  facet_row_q <- sym(facet_row)
  facet_col_q <- sym(facet_col)

  coverage_labeller <- as_labeller(function(x) paste0("Coverage: ", x))

  ggplot(data, aes(x = !!xq, y = !!yq, color = !!colorq, shape = !!shapeq)) +
    geom_jitter(width = 0.2, height = 0, size = 5, alpha = 0.9) +
    facet_grid(rows = vars(!!facet_row_q), cols = vars(!!facet_col_q),
               labeller = labeller(
                 !!facet_col := label_value,
                 !!facet_row := coverage_labeller
               )) +
    scale_shape_manual(values = shape_values) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 16, face = "plain"),
      strip.text.x = element_text(size = 16, face = "bold"),
      strip.text.y = element_text(size = 16, face = "plain")
    ) +
    labs(
      title = title,
      x = x,
      y = y,
      color = color_var,
      shape = shape_var
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))
}



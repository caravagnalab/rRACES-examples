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
    mutate(
      GroundTruth = Signature %in% gt_sigs,
      SparseSignatures = Signature %in% ss_sigs,
      SigProfiler = Signature %in% sp_sigs
    )
  
  long_df <- presence_df %>%
    pivot_longer(cols = c(GroundTruth, SparseSignatures, SigProfiler),
                 names_to = "Method",
                 values_to = "Present") %>%
    filter(Present) %>%
    select(-Present) %>%
    mutate(Method = factor(Method, levels = c("GroundTruth", "SparseSignatures", "SigProfiler")),
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


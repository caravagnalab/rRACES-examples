# Map SparseSignatures to COSMIC #

decompose_sparse_signatures <- function(sparsesig_out, cosmic_signatures_path, mut_counts, threshold = 0.1) {

  cosmic_signatures <- read_delim(cosmic_signatures_path) %>%
    column_to_rownames("Type") %>%
    as.matrix()

  de_novo_sign <- t(sparsesig_out[["beta"]]) %>% as.matrix()
  sparsesig_exposures <- sparsesig_out[["alpha"]] %>% as.matrix()

  # Load sample mutation counts
  sample_catalog <- read.table(mut_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(sample_catalog) <- sample_catalog$MutationType
  sample_catalog <- sample_catalog[, -which(names(sample_catalog) == "MutationType")] %>%
    as.matrix()

  # Decompose de novo signatures to COSMIC
  decompose_cosmic <- MutationalPatterns::fit_to_signatures(de_novo_sign, cosmic_signatures)
  sample_fit <- MutationalPatterns::fit_to_signatures(sample_catalog, cosmic_signatures)

  sample_exposure <- sample_fit$contribution
  cosmic_contribution <- decompose_cosmic$contribution

  # Get all SBS mixtures and significant SBS signatures
  all_sbs_mixtures <- get_all_sbs_mixtures(cosmic_contribution)
  all_significant_sbs <- get_all_significant_sbs(cosmic_contribution, threshold = threshold)

  # Get the unique SBS names and filter the sample exposure
  unique_sbs_names <- get_unique_sbs_names(all_significant_sbs)
  filtered_contribution <- t(sample_exposure[rownames(sample_exposure) %in% unique_sbs_names, ,drop = FALSE])

  return(list(
    filtered_contribution = filtered_contribution,
    cosmic_contribution = cosmic_contribution,
    all_sbs_mixtures = all_sbs_mixtures,
    all_significant_sbs = all_significant_sbs,
    unique_sbs_names = unique_sbs_names
  ))
}

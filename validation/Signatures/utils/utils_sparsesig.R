# Map SparseSignatures to COSMIC #

library("MutationalPatterns")


get_sbs_mixture <- function(cosmic_contribution, de_novo_sig) {
  if (is.numeric(de_novo_sig)) {
    sbs_mixture <- cosmic_contribution[, de_novo_sig, drop = FALSE]
  } else {
    if (!(de_novo_sig %in% colnames(cosmic_contribution))) {
      stop("Specified de novo signature not found in contribution matrix")
    }
    sbs_mixture <- cosmic_contribution[, de_novo_sig, drop = FALSE]
  }
  sbs_mixture <- sbs_mixture[sbs_mixture > 0, , drop = FALSE]

  return(sbs_mixture)
}


get_all_sbs_mixtures <- function(cosmic_contribution) {
  de_novo_sigs <- colnames(cosmic_contribution)

  all_sbs_mixtures <- sapply(de_novo_sigs, function(sig) {
    get_sbs_mixture(cosmic_contribution, sig)
  }, simplify = FALSE)

  return(all_sbs_mixtures)
}


get_significant_sbs <- function(cosmic_contribution, de_novo_sig, threshold) {
  sbs_mixture <- get_sbs_mixture(cosmic_contribution, de_novo_sig)
  significant_sbs <- sbs_mixture[sbs_mixture >= threshold, , drop = FALSE]

  return(significant_sbs)
}


get_all_significant_sbs <- function(cosmic_contribution, threshold) {
  de_novo_sigs <- colnames(cosmic_contribution)

  all_significant_sbs <- lapply(de_novo_sigs, function(sig) {
    get_significant_sbs(cosmic_contribution, sig, threshold)
  })

  names(all_significant_sbs) <- de_novo_sigs

  return(all_significant_sbs)
}


get_unique_sbs_names <- function(all_significant_sbs) {
  unique_sbs <- unique(unlist(lapply(all_significant_sbs, function(data) {
    if (!is.null(data) && nrow(data) > 0) {
      return(rownames(data))
    } else {
      return(NULL)
    }
  })))

  return(unique_sbs)
}



decompose_sparsesig_toCosmic <- function(sparsesig_out, cosmic_signatures_path, mut_counts, threshold = 0.2) {

  cosmic_signatures <- read_delim(cosmic_signatures_path) %>%
    column_to_rownames("Type") %>%
    as.matrix()

  de_novo_sign <- t(sparsesig_out[["beta"]]) %>% as.matrix()
  de_novo_sign <- de_novo_sign %>% as.data.frame() %>% dplyr::select(-Background) %>% as.matrix()
  sparsesig_exposures <- subset(sparsesig_out[["alpha"]], select = -Background) %>% as.matrix()

  # Load sample mutation counts
  sample_catalog <- mut_counts
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

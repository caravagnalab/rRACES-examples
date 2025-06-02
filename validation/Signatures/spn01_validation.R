setwd("/orfeo/cephfs/scratch/cdslab/kdavydzenka/ProCESS-examples/validation/Signatures")

pkgs <- c("tidyverse", "ggplot2", "caret", "ggtext", "reshape2", "lsa", "Metrics", "pheatmap", "MutationalPatterns", "ggvenn")
sapply(pkgs, require, character.only = TRUE)

source("utils/utils_getters.R")
source("utils/utils_validation.R")

# Get Process exposures
main_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
phylo_forest <- getter_process(MAIN_PATH = main_path,
                                SPN_ID = "SPN01",
                                coverage = 50,
                                purity = 0.6,
                                timepoint = NULL,
                                sample_id = NULL,
                                type = "phylo")

process_exposures <- get_exposures(phylo_forest)
process_exposures <- process_exposures[, c("signature", "exposure")]
process_exposures <- process_exposures[grepl("^SBS", process_exposures$signature), ]



# Get tumourevo data
base_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
paths_signatures <- get_signatures(
  PATH = base_path,
  SPN = "SPN01",
  coverage = 50,
  purity = 0.6,
  dataset = "SCOUT",
  context = "SBS"
)

signature_data <- load_signature_data(paths_signatures)

names(signature_data) <- c("Sparsesig", "Sigprofiler_CosmicExposure", "Sigprofiler_CosmicSignature",
                           "Sigprofiler_DenovoExposure", "Sigprofiler_DenovoSignature", "mut_counts")


## Map De novo SparseSignatures to COSMIC using cosine similarity ##

# Extract de novo signatures and exposures
de_novo_signatures <- t(sparsesig_output[["beta"]]) %>% as.matrix()
de_novo_exposures <-  sparsesig_output[["alpha"]]
colnames(de_novo_exposures) <- c("Background", "Denovo_1")

# Load COSMIC reference
cosmic_path <- "COSMIC_v3.4/COSMIC_v3.4_SBS_GRCh38.txt"
cosmic_signatures <- read.delim(cosmic_path) %>%
  column_to_rownames("Type") %>%
  as.matrix()

#Compute cosine similarity
similarity_matrix <- MutationalPatterns::cos_sim_matrix(de_novo_signatures, cosmic_signatures)

# Remap exposures proportional to cosine similarity
similarity_matrix_noBackground <- similarity_matrix[rownames(similarity_matrix) != "Background", , drop = FALSE]

threshold = 0.5

best_matches <- apply(similarity_matrix_noBackground, 1, function(similarities) {
  above_idx <- which(similarities >= threshold)
  if (length(above_idx) == 0) {
    return(NA)  # No matches
  } else {
    matches <- similarities[above_idx]
    matches <- sort(matches, decreasing = TRUE)  
    return(matches)
  }
})


best_matches_denovo <- best_matches %>% as.data.frame()
rownames(best_matches_denovo) <- "Denovo_1"

background_row <- as.data.frame(similarity_matrix["Background", , drop = FALSE])
max_col <- which.max(background_row[1, ])
background_max <- background_row[, max_col, drop = FALSE]

sim_matrix_all <- list(best_matches_denovo, background_max)
names(sim_matrix_all) <- c("Denovo_1", "Background")

cosmic_sigs <- unique(unlist(lapply(sim_matrix_all, names)))
samples <- rownames(de_novo_exposures)
remapped_exposures <- matrix(0, nrow = length(samples), ncol = length(cosmic_sigs),
                             dimnames = list(samples, cosmic_sigs))


for (de_novo_sig in colnames(de_novo_exposures)) {
  if (de_novo_sig == "Background") {
    for (sample in samples) {
      remapped_exposures[sample, "SBS5"] <- remapped_exposures[sample, "SBS5"] +
        de_novo_exposures[sample, de_novo_sig]
    }
  } else {
    sim_weights <- sim_matrix_all[[de_novo_sig]]

    if (!is.null(sim_weights) && length(sim_weights) > 0 && sum(sim_weights) > 0) {
      sim_weights <- as.numeric(sim_weights)
      names(sim_weights) <- names(sim_matrix_all[[de_novo_sig]])  # preserve names
      sim_weights <- sim_weights / sum(sim_weights)

      for (sample in samples) {
        exposure <- de_novo_exposures[sample, de_novo_sig]

        remapped_exposures[sample, names(sim_weights)] <-
          remapped_exposures[sample, names(sim_weights)] +
          exposure * sim_weights
      }
    }
  }
}



# Sum exposures of SBS44 and SBS6 into SBS15
signatures_to_merge <- c("SBS15", "SBS44", "SBS6")
existing_sigs <- intersect(signatures_to_merge, colnames(remapped_exposures))

if (all(c("SBS15") %in% existing_sigs)) {
  remapped_exposures[, "SBS15"] <- rowSums(remapped_exposures[, existing_sigs, drop = FALSE])

  # Remove the merged columns except SBS15
  to_remove <- setdiff(existing_sigs, "SBS15")
  remapped_exposures <- remapped_exposures[, !(colnames(remapped_exposures) %in% to_remove)]
} else {
  warning("SBS15 not found in remapped_exposures columns")
}

# Calculate proportions
row_totals <- rowSums(remapped_exposures)
remapped_exposures_prop <- remapped_exposures / row_totals


## Match each estimated signature to its closest true signature ##
# Sparsesinatures
n_samples <- nrow(remapped_exposures_prop)
sample_names <- rownames(remapped_exposures_prop)
signature_names <- process_exposures$signature
true_exposure_matrix <- matrix(rep(process_exposures$exposure, each = n_samples),
                               nrow = n_samples,
                               dimnames = list(sample_names, signature_names))
true_exposure_df <- as.data.frame(true_exposure_matrix)

true_exposure_df$SBS15 <- rowSums(true_exposure_df[, c("SBS1", "SBS18", "SBS88")])
true_exposure_df_final <- true_exposure_df[, c("SBS5", "SBS15", "SBS17b")]

# Sigprofiler
sigprofiler_exposures <- signature_data[["Sigprofiler_CosmicExposure"]]
sigprofiler_exposures <- sigprofiler_exposures[, !colnames(sigprofiler_exposures) %in% "Samples"]
sigprofiler_exposure_prop <- sigprofiler_exposures / rowSums(sigprofiler_exposures)

n_samples <- nrow(remapped_exposures_prop)
sample_names <- rownames(remapped_exposures_prop)
signature_names <- process_exposures$signature
true_exposure_matrix <- matrix(rep(process_exposures$exposure, each = n_samples),
                               nrow = n_samples,
                               dimnames = list(sample_names, signature_names))
true_exposure_df <- as.data.frame(true_exposure_matrix)

true_exposure_df$SBS15 <- rowSums(true_exposure_df[, c("SBS18", "SBS88")])
true_exposure_df_final <- true_exposure_df[, c("SBS1", "SBS5", "SBS15", "SBS17b")]


## Signatures comparison across methods ##
true_signatures <- colnames(true_exposure_matrix)
sparsesignatures_found <- colnames(remapped_exposures_prop)
sigprofiler_found <- colnames(sigprofiler_exposure_prop)

signature_sets <- list(
  True = true_signatures,
  SparseSignatures = sparsesignatures_found,
  SigProfiler = sigprofiler_found
)

venn_plot <- ggvenn(signature_sets,
       fill_color = c("salmon", "skyblue", "lightgreen"),
       show_percentage = FALSE,
       stroke_size = 0.5,
       set_name_size = 5)


signatures_eval_sparsesig <- evaluate_signatures(true_signatures, sparsesignatures_found)
signatures_eval_sigprofiler <- evaluate_signatures(true_signatures, sigprofiler_found)

conf_matrix <- signatures_eval$ConfusionMatrix
cm_plot_sparsesig <- plot_confusion_matrix(conf_matrix, tool_name = "SparseSignatures")
cm_plot_sigprofiler <- plot_confusion_matrix(conf_matrix, tool_name = "SigProfiler")


## Compute exposure similarities ##
# Sum SBS1 into SBS15

if (all(c("SBS1", "SBS15") %in% colnames(sigprofiler_exposure_prop))) {
  sigprofiler_exposure_prop[, "SBS15"] <- sigprofiler_exposure_prop[, "SBS15"] +
    sigprofiler_exposure_prop[, "SBS1"]

  sigprofiler_exposure_prop <- sigprofiler_exposure_prop[, setdiff(colnames(sigprofiler_exposure_prop), "SBS1")]

} else {
  warning("One of the signatures (SBS1 or SBS15) is missing in the data.")
}



if (all(c("SBS1", "SBS15") %in% colnames(true_exposure_df_final))) {

  true_exposure_df_final[, "SBS15"] <- true_exposure_df_final[, "SBS15"] +
    true_exposure_df_final[, "SBS1"]

  true_exposure_df_final <- true_exposure_df_final[, setdiff(colnames(true_exposure_df_final), "SBS1")]

} else {
  warning("Either SBS1 or SBS15 does not exist in the true exposure data.")
}

results_exposures_val <- validate_exposures(
  process_exposures = true_exposure_df_final,
  sparsesig_out = sparsesig_exp,
  sigprofiler_out = sigprofiler_exp
)


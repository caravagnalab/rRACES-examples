setwd("/orfeo/cephfs/scratch/cdslab/kdavydzenka/ProCESS-examples/validation/Signatures")

pkgs <- c("tidyverse", "ggplot2", "caret", "ggtext", "reshape2", "lsa", "Metrics", "pheatmap", "MutationalPatterns")
sapply(pkgs, require, character.only = TRUE)

source("utils/utils_getters.R")
source("utils/utils_sparsesig.R")
source("utils/utils_validation.R")

# Get Process data
main_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
phylo_forest <- getter_process(MAIN_PATH = main_path,
                                SPN_ID = "SPN03",
                                coverage = 50,
                                purity = 0.6,
                                timepoint = NULL,
                                sample_id = NULL,
                                type = "phylo")

process_exposures <- get_exposures(phylo_forest)
process_exposures <- process_exposures %>% 
  remove_rownames() %>% 
  column_to_rownames(var="sample")


# Get tumourevo data
base_path <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
paths_signatures <- get_signatures(
  PATH = base_path,
  SPN = "SPN03",
  coverage = 50,
  purity = 0.6,
  dataset = "SCOUT",
  context = "SBS"
)

signature_data <- load_signature_data(paths_signatures)

names(signature_data) <- c("Sparsesig", "Sigprofiler_CosmicExposure", "Sigprofiler_CosmicSignature",
                           "Sigprofiler_DenovoExposure", "Sigprofiler_DenovoSignature", "mut_counts")


# Map SparseSignatures to COSMIC #

sparsesig_output <- signature_data[["Sparsesig"]]
cosmic_path <- "COSMIC_ref/COSMIC_v3.4_SBS_GRCh38.txt"
mut_counts <- signature_data[["mut_counts"]]

res_decomposition <- decompose_sparsesig_toCosmic(sparsesig_out = sparsesig_output,
                                      cosmic_signatures_path = cosmic_path,
                                      mut_counts = mut_counts,
                                      threshold = 0.33)

sparsesig_backgroud <- sparsesig_exp %>%
  dplyr::select(Background) %>%
  dplyr::rename(SBS5 = Background)

sparsesig_cosmic_denovo <- res_decomposition[["filtered_contribution"]]
sparsesig_exposures <- cbind(sparsesig_backgroud, sparsesig_cosmic_denovo)
rownames(sparsesig_exposures) <- rownames(process_exposures)



# Confusion matrix

signatures_validation <- generate_confusion_matrix(process_sign = exposures,
                                             sparsesig_sign = sparsesig_exposures,
                                             sigprofiler_sign = signature_data[["Sigprofiler_CosmicExposure"]])


plot_confusion_matrix(signatures_validation$conf_matrix_sparsesig, "SparseSignatures")
plot_confusion_matrix(signatures_validation$conf_matrix_sigprofiler, "Sigprofiler")

sparsesig_out <- sparsesig_exposures

sigprofiler_out <- signature_data[["Sigprofiler_CosmicExposure"]] %>%
  remove_rownames() %>%
  column_to_rownames(var="Samples")
rownames(sigprofiler_out) <- rownames(process_exposures)


# Exposures validation

results_exposures_val <- validate_exposures(
  process_exposures = process_exposures,
  sparsesig_out = sparsesig_out,
  sigprofiler_out = sigprofiler_out)

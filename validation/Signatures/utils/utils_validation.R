# Compute confusion matrix #

validate_signatures <- function(process_sign, sparsesig_sign, sigprofiler_sign) {
  ground_truth <- colnames(process_sign)[-1]
  sparsesig <- colnames(sparsesig_sign)
  sigprofiler <- colnames(sigprofiler_sign)[-1]
  all_signatures <- unique(c(ground_truth, sparsesig, sigprofiler))
  
  binary_matrix <- function(target, ref) as.numeric(target %in% ref)
  conf_matrix_sparsesig <- table(Predicted = binary_matrix(all_signatures, sparsesig), True = binary_matrix(all_signatures, ground_truth))
  conf_matrix_sigprofiler <- table(Predicted = binary_matrix(all_signatures, sigprofiler), True = binary_matrix(all_signatures, ground_truth))
  
  return(list(conf_matrix_sparsesig = conf_matrix_sparsesig, conf_matrix_sigprofiler = conf_matrix_sigprofiler))
}

library(rRACES)
library(CNAqc)
library(dplyr)
library(ggplot2)

relative_to_absolute_coords_pos = function(x, ref = "GRCh38") {
    reference_genome = CNAqc:::get_reference(ref)
    vfrom = reference_genome$from
    names(vfrom) = reference_genome$chr
    
    # colnames(x) = c("chr", "pos")
    x = x %>%
      dplyr::rename(pos=from) %>% 
      dplyr::mutate(chr = paste0("chr", chr)) %>%
      dplyr::mutate(pos = pos + vfrom[chr]) 
    return(x)
}

plot_DR_n <- function(
    seq_res,
    sample,
    chromosomes = NULL,
    N = 5000) {
  
  data <- rRACES:::get_seq_data(seq_res, sample, chromosomes)
  
  Ntotal <- nrow(data$tumour)
  N = min(N, Ntotal)
  mean_DP_normal <- mean(data$normal$DP)
  mean_DP_tumor <- mean(data$tumour$DP)
  d <- data$tumour %>%
    dplyr::left_join(data$normal, suffix = c(".tumour", ".normal"),
                     by = c("chr", "from", "ref", "alt")) %>%
    dplyr::mutate(DR = (DP.tumour*mean_DP_normal) / (DP.normal*mean_DP_tumor)) %>%
    dplyr::sample_n(N) %>%
    dplyr::arrange(chr, from)
  d_mod <- relative_to_absolute_coords_pos(d,ref = "GRCh38")
  CNAqc:::blank_genome() +
    geom_point(data=d_mod, aes(x = pos, y = DR), size = 0.5) +
    ylim(c(0, 4)) +
    ylab("DR")
}


plot_BAF_n <- function(
    seq_res,
    sample,
    chromosomes = NULL,
    cuts = c(0, 1),
    N = 5000) {
  data <- rRACES:::get_seq_data(seq_res, sample, chromosomes)

  Ntotal <- nrow(data$tumour)
  N = min(N, Ntotal)
  
  d <- data$tumour %>%
   #dplyr::group_by(chr) %>%
   dplyr::sample_n(N) %>%
   #dplyr::ungroup() %>%
   dplyr::arrange(chr, from) %>%
   dplyr::mutate(abs_pos = 1:dplyr::n()) %>%
   dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))

  d_mod <- relative_to_absolute_coords_pos(d,ref = "GRCh38")
  CNAqc:::blank_genome() +
    geom_point(data=d_mod, aes(x = pos, y = VAF), size = 0.5)+
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::labs(x = "", y = "BAF")
}

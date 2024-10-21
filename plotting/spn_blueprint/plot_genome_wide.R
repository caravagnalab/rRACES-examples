library(rRACES)
library(CNAqc)
library(dplyr)
library(ggplot2)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/plotting/spn_blueprint/utils.R")

#' Plot Depth Ratio (DR) Genome-wide normalize
#'
#' This function normalize sequencing depth of tumour and normal samples
#' generates a plot showing the DEpth Ratio (DR) genome-wide for a specific sample.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot. Default is all autosomes and sex chromosomes.
#' @param n The number of mutations to sample for plotting. Default is 5000.
#' @return A ggplot2 object showing the DR distribution across the genome.
#'
#' @export

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

#' Plot B-allele frequency (BAF) Genome-wide normalize
#'
#' This function generates a plot showing the B-allele frequency (DR) genome-wide for a specific sample.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot. Default is all autosomes and sex chromosomes.
#' @param n The number of mutations to sample for plotting. Default is 5000.
#' @return A ggplot2 object showing the DR distribution across the genome.
#'
#' @export

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


genome_wide_plots <- function(x,seq_results,sample_id){
  g_seq <- seq_results %>% filter(classes=="germinal")
  dr <- plot_DR_n(seq_res = g_seq,sample = sample_id)
  baf <- plot_BAF_n(seq_res = g_seq,sample = sample_id)
  vaf_histo <- plot_data_histogram(x, which = 'VAF')
  dp_histo <- plot_data_histogram(x, which = 'DP')
  cov=100
  x_original <- x
  ccf_lables <- unique(x$cna$ccf_label) ### get the different ccf
  print(ccf_lables)
  segs <- list()
  for (c in ccf_lables){
    x$cna <- x$cna %>% filter(ccf_label==c)
    ccf_ratio <- median(x$cna$ratio)
    seg <- plot_segments(x = x,highlight = c("1:0","1:1","2:1","2:2"))+
	    ggplot2::ggtitle(paste0("CCF: ",ccf_ratio))
    segs[[c]] <- seg
    x<-x_original
  }
  p1 <- vaf_histo/dp_histo+plot_layout(design="AAAAA\nAAAAA\nAAAAA\nBBBBB\nBBBBB\nBBBBB")+
	  plot_annotation(title=sample_id,subtitle=
			  paste0("Simulated coverage: ",cov,"\nSimulated purity: ",x$purity))
  p2 <- dr/baf+plot_layout(design="AAAAA\nAAAAA\nBBBBB\nBBBBB")+
	  plot_annotation(title=sample_id,subtitle=
			  paste0("Simulated coverage: ",cov,"\nSimulated purity: ",x$purity))
  p3 <- wrap_plots(segs)+plot_layout(design="AAAAA\nAAAAA\nBBBBB\nBBBBB\nCCCCC\nCCCCC")
  return(list("histo"=p1,"genome_wide"=p2,"segments"=p3))
}

library(rRACES)
library(CNAqc)
library(dplyr)
library(ggplot2)
library(ggalluvial)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/REPO_UPDATED/rRACES-examples/report/plotting/utils.R")

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
  mean_DP_normal <- mean(na.omit(data$normal$DP))
  mean_DP_tumor <- mean(na.omit(data$tumour$DP))
  d <- data$tumour %>%
    dplyr::left_join(data$normal, suffix = c(".tumour", ".normal"),
                     by = c("chr", "from", "ref", "alt")) %>%
    dplyr::mutate(DR = (DP.tumour*mean_DP_normal) / (DP.normal*mean_DP_tumor)) %>%
    dplyr::sample_n(N) %>%
    dplyr::arrange(chr, from) %>% 
    dplyr::filter(causes.tumour=="" & causes.normal=="")
  d_mod <- relative_to_absolute_coords_pos(d,ref = "GRCh38")
  CNAqc:::blank_genome() +
    geom_point(data=d_mod, aes(x = pos, y = DR), size = 0.5) +
    ylim(c(0, 4)) +
    ylab("DR")+
    ggplot2::ggtitle(sample)
}

#' Plot B-allele frequency (BAF) Genome-wide normalize
#'
#' This function generates a plot showing the B-allele frequency (BAF) genome-wide for a specific sample.
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
  data <- rRACES:::get_seq_data(as.data.frame(seq_res), sample, chromosomes)

  Ntotal <- nrow(data$tumour)
  N = min(N, Ntotal)
  # cov = mean(na.omit(data$tumour$DP))
  d <- data$tumour %>%
   #dplyr::group_by(chr) %>%
   dplyr::sample_n(N) %>%
   #dplyr::ungroup() %>%
   dplyr::arrange(chr, from) %>%
   dplyr::mutate(abs_pos = 1:dplyr::n()) %>%
   dplyr::filter(VAF <= max(cuts), VAF >= min(cuts)) %>% 
  dplyr::filter(causes!=" + errors")
  color_values <- c("50" = "lightgrey", "100" = "darkgrey", "150" = "grey", "200" = "black")
  d_mod <- relative_to_absolute_coords_pos(d,ref = "GRCh38")
  CNAqc:::blank_genome() +
    geom_point(data=d_mod, aes(x = pos, y = VAF), size = 0.5)+
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::labs(x = "", y = "BAF")+
    # scale_fill_manual(values = color_values[as.character()]
    ggplot2::ggtitle(sample)
}


#' Plot variant allele frequency (VAF) Genome-wide normalize
#'
#' This function generates a plot showing the VAF genome-wide for a specific sample.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot. Default is all autosomes and sex chromosomes.
#' @param n The number of mutations to sample for plotting. Default is 5000.
#' @return A ggplot2 object showing the DR distribution across the genome.
#'
#' @export

plot_VAF_n <- function(
    seq_res,
    sample,
    chromosomes = NULL,
    cuts = c(0.01, 1),
    N = 80000) {
  data <- rRACES:::get_seq_data(as.data.frame(seq_res), sample, chromosomes)
  
  Ntotal <- nrow(data$tumour)
  N = min(N, Ntotal)
  # cov = mean(na.omit(data$tumour$DP))
  d <- data$tumour %>%
    #dplyr::group_by(chr) %>%
    dplyr::sample_n(N) %>%
    #dplyr::ungroup() %>%
    dplyr::arrange(chr, from) %>%
    dplyr::mutate(abs_pos = 1:dplyr::n()) %>%
    dplyr::filter(VAF != 0) %>% 
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts)) %>% 
    dplyr::filter(causes!=" + errors")
  color_values <- c("50" = "lightgrey", "100" = "darkgrey", "150" = "grey", "200" = "black")
  d_mod <- relative_to_absolute_coords_pos(d,ref = "GRCh38")
  CNAqc:::blank_genome() +
    geom_point(data=d_mod, aes(x = pos, y = VAF), size = 0.5)+
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::labs(x = "", y = "BAF")+
    # scale_fill_manual(values = color_values[as.character()]
    ggplot2::ggtitle(sample)
}

plot_VAF_data_n = function(x,
                         karyotypes = c("1:0", '1:1', '2:0', '2:1', '2:2'),
			 cuts)
{
  stopifnot(inherits(x, 'cnaqc'))

  # VAF
  # raw_muts = x$mutations %>%
  #   dplyr::mutate(karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"))

  raw_muts = Mutations(x) %>%
    dplyr::mutate(
      karyotype = ifelse(karyotype %in% karyotypes, karyotype, "other"),
      cna = paste("CNA", cna)
      ) %>%
    dplyr::filter(VAF <= max(cuts), VAF > min(cuts))
 
  colors = CNAqc:::get_karyotypes_colors(unique(raw_muts$karyotype))
  colors['subclonal'] = ggplot2::alpha('purple4', .7)

  ggplot2::ggplot(data = raw_muts, ggplot2::aes(VAF, fill = karyotype)) +
    ggplot2::geom_histogram(binwidth = 0.01) +
    ggplot2::xlim(-0.01, 1.01) +
    CNAqc:::my_ggplot_theme() +
    ggplot2::labs(title = "VAF",
                  caption = paste0(
                    "n = ",
                    nrow(raw_muts),
                    "; VAF < ",min(cuts)," = ",
                    sum(x$mutations$VAF < min(cuts))
                  )) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::facet_wrap(type ~ cna, scales = 'free_y')
}


#' Plot genome wide plots for each sample.
#'
#' This function will wrap-up all the plots at the genome wide level for each sample.
#' The included plots are: Normalized depth ratio, B-allele frequency, VAF and DP
#' distributions for somatic mutations, and allelic copy number events for different
#' identified subclones.
#'
#' @param x A CNAqc object
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample_id  The name of the sample for which the plot will be generated.
#' @param coverage Sequencing depth of the simulation.
#'
#' @export

genome_wide_plots <- function(x,seq_results,sample_id){
  g_seq <- seq_results %>% filter(classes=="germinal")
  dr <- plot_DR_n(seq_res = g_seq,sample = sample_id)
  baf <- plot_BAF_n(seq_res = g_seq,sample = sample_id)
  vaf_histo <- plot_VAF_data_n(x,cuts=c(0.02,1))
  dp_histo <- plot_data_histogram(x, which = 'DP')
  vaf_histo_karyotype <-plot_VAF_data_n(x,cuts=c(0.02,1))+
                facet_wrap(~ karyotype, ncol=3,scales = "free")+
		ggplot2::ggtitle("Variant Allele Spectra per karyotype")
  x_original <- x
  ccf_lables <- unique(x$cna$ccf_label) ### get the different ccf
  print(ccf_lables)
  segs <- list()
  for (c in ccf_lables){
    x$cna <- x$cna %>% filter(ccf_label==c)
    ccf_ratio <- round(median(x$cna$ratio),2)
    seg <- plot_segments(x = x,highlight = c("1:0","1:1","2:1","2:2"))+
	    ggplot2::ggtitle(paste0("Median fraction of tumour cells: ",ccf_ratio))
    segs[[c]] <- seg
    x<-x_original
  }
  p1 <- vaf_histo/dp_histo+plot_layout(design="AAAAA\nAAAAA\nAAAAA\nBBBBB\nBBBBB\nBBBBB")
  p2 <- dr/baf+plot_layout(design="AAAAA\nAAAAA\nBBBBB\nBBBBB")
  if (length(segs)<3){
    n<-length(segs)+1
    for (i in c(n:3)){
       segs[[i]]<-ggplot()
    }
  }
  p3 <- (wrap_plots(segs)+vaf_histo_karyotype)+
	  plot_layout(design="AAAAA\nAAAAA\nBBBBB\nBBBBB\nCCCCC\nCCCCC\nDDDDD\nDDDDD\nDDDDD\nDDDDD")
  return(list("histograms"=p1,"genome_wide"=p2,"segments"=p3))
}

plot_clone_segments <- function(files_cna){
  offset <- 0.05
  sample_names <- sapply(files_cna, function(path) {
    base_name <- basename(path)
    sub("_cna.rds$", "", base_name)
  })

  upper <- 0.89
  lower <- 0.11
  
  cna_seg <- lapply(files_cna, function(f){
    base_name <- basename(f)
    sample <- sub("_cna.rds$", "", base_name)
    readRDS(f) %>% 
      mutate(CN_type = ifelse(ratio < upper & ratio > lower, 'sub-clonal', 'clonal'),
             CN = paste(major, minor, sep = ':'),
             seg_id = paste(chr,begin,end, sep = ':'),
             sample = sample)
  })
  names(cna_seg) <- sample_names
  
  data <- absolute_to_relative_coordinates(cna_seg %>% bind_rows() %>% mutate(chr = paste0('chr', chr))) %>% 
    filter(!(CN_type == 'sub-clonal' & ratio  > upper)) %>% 
    filter(!(CN_type == 'clonal' & ratio <= lower))  %>% 
    mutate(ratio = ifelse(ratio < upper & ratio > lower, ratio, 1)) %>% 
    mutate(ratio = round(ratio, digits=1))
  
  plt <- blank_genome() +
    geom_rect(data = data, aes(xmin = begin, xmax = end, ymin = -Inf, ymax = Inf, fill = factor(CN)), alpha = 0.3) +
    geom_segment(data = data, aes(x = begin, xend = end, y = major+offset, yend = major+offset), col = 'red', size = 1) +
    geom_segment(data = data, aes(x = begin, xend = end, y = minor-offset, yend = minor-offset), col = 'steelblue', size = 1) +
    scale_fill_manual(values = get_karyotypes_colors(unique(data$CN))) + 
    ylab('CN') + 
    xlab('position') + 
    ggplot2::guides(fill = ggplot2::guide_legend('CN', override.aes = list(alpha = 1))) + 
    ylim(-0.2, 4.2) + 
    ggh4x::facet_nested(sample + CN_type + ratio ~.)
  return(plt)
}


alluvial_plot_karyotypes <- function(files_cna){
  
  samples <- sapply(files_cna, function(path) {
    base_name <- basename(path)
    sub("_cna.rds$", "", base_name)
  })
  bulk_cell_cna_list <- lapply(files_cna, function(x){
    base_name <- basename(x)
    sample <- sub("_cna.rds$", "", base_name)
    c <- readRDS(x) %>% 
      mutate(sample=sample)
  })
  
  cna_data <- do.call("rbind",bulk_cell_cna_list)
  karyotype_of_interest <- c("1:0","1:1","2:1","2:0","2:2","3:1","3:2")
  
  cna_data_summary <- cna_data %>% 
    filter(ratio >=0.1) %>% 
    mutate(karyotype=paste(major,minor,sep=":")) %>%
    mutate(karyotype=case_when(!(karyotype %in%karyotype_of_interest)~"Others",
                               karyotype %in%karyotype_of_interest ~ karyotype)) %>% 
    group_by(sample,karyotype) %>% 
    summarise(n = n()) 
  
  
  all_combinations <- expand.grid(
    sample = unique(cna_data_summary$sample),
    # mutant = unique(cna_data_summary$mutant),
    karyotype = unique(cna_data_summary$karyotype),
    stringsAsFactors = FALSE
  )
  
  df_full <- all_combinations %>%
    left_join(cna_data_summary, by = c("sample","karyotype")) %>%
    mutate(n = ifelse(is.na(n), 0, n))  # Replace NA values with 0 for missing combinations
  
  df_fractions <- df_full %>%
    group_by(sample) %>%
    mutate(total_n = sum(n),
           fraction = n / total_n*100) %>%
    select(-total_n)                
  
  plot <-df_fractions %>% 
    mutate(clone=rep(c("clone_1","clone_2"),length.out = n())) %>% 
    ggplot(aes(x = sample, stratum = karyotype, alluvium = karyotype, y = fraction, fill=karyotype)) +
    geom_alluvium(aes(fill = karyotype), alpha = 0.5) +  # Flows of karyotypes
    geom_stratum() +  # Create clone blocks in each sample
    # geom_text(stat = "stratum", aes(label = clone), size = 5, color = "white") +  # Clone labels
    my_ggplot_theme() +
    scale_fill_manual(values=get_karyotypes_colors(unique(df_full$karyotype)))+
    labs(x = "Sample", y = "Fraction of Segments") +
    scale_y_continuous(expand = c(0, 0))  # Keep Y-axis clean
  return(plot)
}

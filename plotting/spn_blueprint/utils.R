rm(list = ls())
library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)
#library(CNAqc)
#source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/scripts/my_functions/plot_genome_wide.R")
#'
#' Convert realtive coordinates to absolute coordinates
#'
#' This function takes CNAqc bject and converts relative coordinates
#' to absolute coordinates given an input reference genome.
#'
#' @param x CNAqc object
#' @param ref name of the reference genome
#' @return x CNAqc object with absoluute coordinates
#'
#' @export


relative_to_absolute_coords_pos = function(x, ref = "GRCh38") {
  reference_genome = CNAqc:::get_reference(ref)
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  x = x %>%
      dplyr::rename(pos=from) %>%
      dplyr::mutate(chr = paste0("chr", chr)) %>%
      dplyr::mutate(pos = pos + vfrom[chr])
  return(x)
}


signatures_palette <- function(phylo_forest,seed){
  ref_path <- phylo_forest$get_reference_path()
  SBS_table_path <- gsub(pattern = "reference.fasta",replacement = "SBS_signatures.txt",x = ref_path)
  IDS_table_path <- gsub(pattern = "reference.fasta",replacement = "indel_signatures.txt",x = ref_path)
  SBS_table <-  read.csv(SBS_table_path,header=T,sep="\t")
  SBS_sign <- colnames(SBS_table)
  IDS_table <-  read.csv(IDS_table_path,header=T,sep="\t")
  IDS_sign <- colnames(IDS_table)
  set.seed(seed)
  return(Polychrome::createPalette(length(sigs), c("#6B8E23","#4169E1"), M=1000,
                                   target="normal", range=c(15,80)) %>%
           setNames(sigs))
}
#'
#' Generates CNAqc object given simulate_seq() dataframe and
#' phylogenetic forest object. The final result will be a CNAqc 
#' object.
#'

races2cnaqc <- function(seq_results,phylo_forest,sample_id,ref,purity){
  ref_path <- phylo_forest$get_reference_path()
  driver_table_path <- gsub(pattern = "reference.fasta",replacement = "drivers.txt",x = ref_path)
  driver_table <-  read.csv(driver_table_path,header=T,sep="\t")
  known_drivers <- driver_table %>%
    dplyr::mutate(chr=as.character(chr)) %>%
    dplyr::rename(driver_label=driver_code)
  bulk <- phylo_forest$get_bulk_allelic_fragmentation(sample_id)
  cna <- bulk %>% dplyr::rename("Major"=major,"from"=begin,"to"=end) %>%
    dplyr::filter(ratio>=0.08)
  cna <- cna %>%
    group_by(chr, from, to) %>%            # Group by chr, begin, and end
    arrange(desc(ratio)) %>%                 # Arrange by 'ratio' in descending order
    mutate(ccf_label = paste0("ccf_", rank(-ratio)))  # Assign rank based on 'ratio'

  mutations <- rRACES::seq_to_long(seq_results) %>%
    dplyr::filter(sample_name==sample_id & classes!="germinal") %>%
    dplyr::filter(VAF!=0) %>% mutate(is_driver=FALSE) %>%
    left_join(known_drivers,by=c("chr","from","to","ref","alt")) %>%
    dplyr::mutate(
      is_driver = ifelse(!is.na(driver_label), TRUE,
                         ifelse(is.na(driver_label) & classes == "driver", TRUE, FALSE)),
      driver_label = ifelse(is.na(driver_label) & classes == "driver",
                            paste(chr, from, ref, alt, sep = ":"), driver_label))
  x <- CNAqc::init(mutations = mutations,cna = cna,
                   purity = purity,sample = sample_id,
                   ref = ref)
  return(x)
}


get_classes_colors <- function(classes){
  color = c(`driver` = "firebrick4",`passenger` = ggplot2::alpha("tan2",0.4),`pre-neoplastic` = "cornflowerblue",
            `germinal` = "darkolivegreen")
  missing = setdiff(names(color), classes)
  nmissing = length(missing)
  c(color, CNAqc:::nmfy(missing, rep("gray", nmissing)))
}



get_legend <- function(col_palette){
  df <- data.frame(type = names(col_palette), color = col_palette)
  p <- ggplot(df, aes(x = type, fill = type)) +
    geom_bar() +
    scale_fill_manual(values = col_palette) +
    theme_void() +  # Remove axes and background
    guides(fill = guide_legend(title = "Classes & Causes"))

  legend_plot <- ggpubr::get_legend(p,position = "right")
  pl <- ggpubr::as_ggplot(legend_plot)
  return(pl)
}


squareplot = function(seq_res, samples_list,chrom)
{
  row_plots = NULL
  for (s in seq(samples_list))
  {
    sn = samples_list[s]
    s_seq <- seq_res %>% filter(classes!="germinal")
    s_seq_long <- s_seq %>% rRACES::seq_to_long()
    plot_vaf <- s_seq_long %>%
      filter(sample_name==sn & chr==chrom) %>%
      filter(VAF!=0) %>%
      ggplot(aes(x=VAF)) +geom_histogram(binwidth = 0.01) +
      xlim(c(0,1))+
      ggplot2::ggtitle(label = sn) +
      CNAqc:::my_ggplot_theme()



    mb = list(plot_vaf+ labs(title = sn) )

    idx_pre = 1:s
    idx_post = s:length(samples_list)

    pl_r = pl_l = NULL
     
    #palette <- RColorBrewer::brewer.pal(n = length(unique(s_seq$causes)), name = "Set3")
    #col_causes <- setNames(palette, unique(s_seq$causes))
    #cols_causes <- rRACES:::get_colors_for(unique(s_seq$causes))
    #col_classes <- c("passenger" = "#CCCC99",
    #                 "pre-neoplastic" = "#006699",
    #                 "driver" = "#990033")
    #cols <- c(col_causes,col_classes)

    if (length(idx_pre) > 1)
      pl_r = lapply(setdiff(idx_pre, s), function(x) {
        s_sn <- s_seq_long %>% filter(sample_name==sn & chr==chrom)
        s_sn_x <- s_seq_long %>% filter(sample_name==samples_list[x] & chr==chrom)
        joined <- full_join(s_sn,s_sn_x,by=c("chr","from","ref","alt","to","causes","classes"))
        plot <- joined %>% ggplot(aes(x=VAF.x,y=VAF.y,col=classes)) + geom_point() +
          CNAqc:::my_ggplot_theme()
          #scale_color_manual(values = col_classes)
        plot + ggplot2::geom_point(alpha = 0.7) +
          ggplot2::xlim(c(-0.01, 1.01)) +
          ggplot2::ylim(c(-0.01, 1.01)) +
          ggplot2::labs(x = sn, y = samples_list[x])+
          ggplot2::theme(legend.position = "none")
      })

    if (length(idx_post) > 1)
      pl_l = lapply(setdiff(idx_post, s), function(x) {
        s_sn <- s_seq_long %>% filter(sample_name==sn & chr==chrom)
        s_sn_x <- s_seq_long %>% filter(sample_name==samples_list[x] & chr==chrom)
        joined <- full_join(s_sn,s_sn_x,by=c("chr","from","ref","alt","to","causes","classes"))
        plot <- joined %>% ggplot(aes(x=VAF.x,y=VAF.y,col=causes)) + geom_point() +
          CNAqc:::my_ggplot_theme()
          #scale_color_manual(values = col_causes)
        plot + ggplot2::geom_point(alpha = 0.7) +
          ggplot2::xlim(c(-0.01, 1.01)) +
          ggplot2::ylim(c(-0.01, 1.01)) +
          ggplot2::labs(x = sn, y = samples_list[x])+
          ggplot2::theme(legend.position = "none")
      })

    plotlist = append(append(pl_r, mb), pl_l)
    row_plot = patchwork::wrap_plots(plotlist)+
      patchwork::plot_layout(guides = "collect",ncol = length(pl_r) + length(pl_l) + 1,nrow = 1)
    row_plots = append(row_plots, list(row_plot))
  }
  #pl <- get_legend(cols)
  patchwork::wrap_plots(row_plots)+
    patchwork::plot_layout(design = "AAAA\nBBBB\nCCCC")+
    patchwork::plot_annotation(title = paste0("Chromosome ", chrom))
}

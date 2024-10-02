rm(list = ls())
library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CNAqc)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/scripts/my_functions/plot_DR_norm.R")

races2cnaqc <- function(seq_results,phylo_forest,sample_id,ref,purity){
  ref_path <- phylo_forest$get_reference_path()
  driver_table_path <- gsub(pattern = "reference.fasta",replacement = "drivers.txt",x = ref_path)
  driver_table <-  read.csv(driver_table_path,header=T,sep="\t")
  known_drivers <- driver_table %>%
    dplyr::mutate(chr=as.character(chr)) %>%
    dplyr::rename(driver_label=driver_code)
  bulk <- phylo_forest$get_bulk_allelic_fragmentation(sample_id)
  cna <- bulk %>% dplyr::rename("Major"=major,"from"=begin,"to"=end) %>%
    dplyr::filter(ratio>=0.05)
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

genome_wide_plots <- function(x,seq_results,sample_id){
  g_seq <- seq_results %>% filter(classes=="germinal")
  seg <- plot_segments(x = x,highlight = c("1:0","1:1","2:1","2:2"))
  dr <- plot_DR_n(seq_res = g_seq,sample = sample_id)
  baf <- plot_BAF_n(seq_res = g_seq,sample = sample_id)
  histo <- ggpubr::ggarrange(
    plot_data_histogram(x, which = 'VAF'),
    plot_data_histogram(x, which = 'DP'),
    plot_data_histogram(x, which = 'NV'),
    ncol = 3,
    nrow = 1
  )
  cov=100
  p = histo/dr/baf/seg+plot_layout(design="AAAA\nAAAA\nBBBB\nCCCC\nDDDD")+
	  plot_annotation(title=sample_id,subtitle=
			  paste0("Simulated coverage: ",cov,"\nSimulated purity: ",x$purity))
  p
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
    plot_vaf <- s_seq %>% rRACES::seq_to_long() %>%
      filter(sample_name==sn & chr==chrom) %>%
      filter(VAF!=0) %>%
      ggplot(aes(x=VAF)) +geom_histogram(binwidth = 0.01) +
      ggplot2::ggtitle(label = sn) +
      CNAqc:::my_ggplot_theme()



    mb = list(plot_vaf+ labs(title = sn) )

    idx_pre = 1:s
    idx_post = s:length(samples_list)

    pl_r = pl_l = NULL
    
    palette <- RColorBrewer::brewer.pal(n = lenght(unique(s_seq$causes)), name = "Set3")
    cols <- setNames(palette, unique(s_seq$causes))
    #cols <- rRACES:::get_colors_for(unique(s_seq$causes))
    col_classes <- c("passenger" = "#CCCC99",
                     "pre-neoplastic" = "#006699",
                     "driver" = "#990033")
    cols <- c(cols,col_classes)

    if (length(idx_pre) > 1)
      pl_r = lapply(setdiff(idx_pre, s), function(x) {
        s_sn <- s_seq %>% rRACES::seq_to_long() %>% filter(sample_name==sn & chr==chrom)
        s_sn_x <- s_seq %>% rRACES::seq_to_long() %>% filter(sample_name==samples_list[x] & chr==chrom)
        joined <- full_join(s_sn,s_sn_x,by=c("chr","from","ref","alt","to","causes","classes"))
        plot <- joined %>% ggplot(aes(x=VAF.x,y=VAF.y,col=classes)) + geom_point() +
          CNAqc:::my_ggplot_theme()+
          scale_color_manual(values = col_classes)
        plot + ggplot2::geom_point(alpha = 0.7) +
          ggplot2::xlim(c(-0.01, 1.01)) +
          ggplot2::ylim(c(-0.01, 1.01)) +
          ggplot2::labs(x = sn, y = samples_list[x])+
          ggplot2::theme(legend.position = "none")
          # ggplot2::theme_bw() +
          # ggplot2::theme(legend.position = "bottom")
        # VIBER::plot_2D(viber_fit_bottom, d1 = sn, d2 = samples_list[x], colors = colors_bottom)
      })

    if (length(idx_post) > 1)
      pl_l = lapply(setdiff(idx_post, s), function(x) {
        s_sn <- s_seq %>% rRACES::seq_to_long() %>% filter(sample_name==sn & chr==chrom)
        s_sn_x <- s_seq %>% rRACES::seq_to_long() %>% filter(sample_name==samples_list[x] & chr==chrom)
        joined <- full_join(s_sn,s_sn_x,by=c("chr","from","ref","alt","to","causes","classes"))
        plot <- joined %>% ggplot(aes(x=VAF.x,y=VAF.y,col=causes)) + geom_point() +
          CNAqc:::my_ggplot_theme() +
          scale_color_manual(values = col_signatures)
        plot + ggplot2::geom_point(alpha = 0.7) +
          ggplot2::xlim(c(-0.01, 1.01)) +
          ggplot2::ylim(c(-0.01, 1.01)) +
          ggplot2::labs(x = sn, y = samples_list[x])+
          ggplot2::theme(legend.position = "none")
          # ggplot2::theme_bw() +
          # ggplot2::theme(legend.position = "bottom")
      })

    # row_plot = cowplot::plot_grid(
    #   plotlist = append(append(pl_r, mb), pl_l),
    #   nrow = 1,
    #   ncol = length(pl_r) + length(pl_l) + 1,
    #   align = 'h',
    #   axis = 'x'
    # )
    plotlist = append(append(pl_r, mb), pl_l)
    row_plot = patchwork::wrap_plots(plotlist)+
      patchwork::plot_layout(guides = "collect",ncol = length(pl_r) + length(pl_l) + 1,nrow = 1)
    row_plots = append(row_plots, list(row_plot))
  }
  pl <- get_legend(cols)
  patchwork::wrap_plots(row_plots)+pl+
    patchwork::plot_layout(design = "AAAD\nBBBD\nCCCD")+
    ggplot2::ggtitle(label=paste0("chr",chrom))

  # cowplot::plot_grid(
  #   plotlist = row_plots,
  #   ncol = 1,
  #   nrow = length(row_plots),
  #   align = 'v'
  # )
}

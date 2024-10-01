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

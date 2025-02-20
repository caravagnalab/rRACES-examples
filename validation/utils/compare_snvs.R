rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
library(vcfR)
library(naniar)
#library(ggExtra)
source("/u/cdslab/ggandolfi/scratch/prj_races/rRACES-examples/races_validation/utils/vcf_parser.R") ## vcf parser

## set ggplot2 theme
my_ggplot_theme<- function (cex = 1) {
    cex_opt = 1
    theme_light(base_size = 10 * cex_opt) + theme(legend.position = "bottom",
        legend.key.size = unit(0.3 * cex_opt, "cm"), panel.background = element_rect(fill = "white"))
}

mutation_calling_accuracy <- function(seq_races,vcf_file,caller, sample_id,top_filter,sample_id_vcf){
  seq_res <- readRDS(seq_races)
  if ("sample_name"%in%colnames(seq_res)){
    seq_res_long<-seq_res
  } else {
	seq_res_long <- seq_to_long(seq_res)
  }

  ## extract sample for which variant calling has been done
  #sample_id <- "SPN01_Sample_1"
  seq_res_long <- seq_res_long %>%
    dplyr::mutate(mutationID=paste0("chr",chr,":",from,":",ref,":",alt)) %>%
    dplyr::filter(classes!="germinal") %>%
    dplyr::filter(sample_name==sample_id) %>%
    dplyr::filter(NV!=0)
  chromsomes <- paste0("chr",unique(seq_res_long$chr))
  cols_to_keep <- c("chr", "from", "to", "ref","alt","mutationID","NV","DP","VAF","FILTER")
  plot_vaf_races <- seq_res_long %>% ggplot(aes(x=VAF)) + geom_histogram(binwidth=0.01)
	  #facet_wrap(~chr,scale="free")
  ## READ vcf
  vcf <- vcfR::read.vcfR(file = vcf_file)
  if (caller=="mutect2"){
    muts <- parse_Mutect(vcf = vcf,filter_mutations = FALSE)
    muts[[sample_id_vcf]]$mutations <-muts[[sample_id_vcf]]$mutations %>%
      dplyr::filter(chr%in%chromsomes) %>%
      dplyr::mutate(mutationID=paste0(chr,":",from,":",ref,":",alt))
    print(muts[[sample_id_vcf]]$mutations[,c(cols_to_keep)])
    merged_df <- merge(seq_res_long, muts[[sample_id_vcf]]$mutations[,c(cols_to_keep)],
                       by = "mutationID", all = TRUE, suffix=c(".races",".caller"))
    plot_vaf_caller <- muts[[sample_id_vcf]]$mutations %>% ggplot(aes(x=VAF)) + geom_histogram(binwidth=0.01)

  } else if (caller=="strelka"){
    muts <- parse_Strelka(vcf = vcf,filter_mutations = FALSE)
    muts$TUMOR$mutations <-muts$TUMOR$mutations %>%
      dplyr::filter(chr%in%chromsomes) %>%
      dplyr::mutate(mutationID=paste0(chr,":",from,":",ref,":",alt))
    merged_df <- merge(seq_res_long, muts$TUMOR$mutations[,c(cols_to_keep)],
                       by = "mutationID", all = TRUE, suffix=c(".races",".caller"))
    plot_vaf_caller <- muts$TUMOR$mutations  %>% ggplot(aes(x=VAF)) + geom_histogram(binwidth=0.01)

  }
  merged_df$type <- ifelse(!is.na(merged_df$VAF.races) & !is.na(merged_df$VAF.caller), "TP",
                           ifelse(is.na(merged_df$VAF.caller), "FN", "FP"))

  ### STATS ###

  matrix_scores <- merged_df %>% count(type)
  tp <- matrix_scores[matrix_scores$type=="TP","n"]
  fn <- matrix_scores[matrix_scores$type=="FN","n"]
  fp <- matrix_scores[matrix_scores$type=="FP","n"]
  sensitivity <- tp/(tp+fn)
  precision <- tp/(tp+fp)
  tp_data <-  merged_df %>% dplyr::filter(type=="TP")
  acc_vaf <- sqrt(mean((tp_data$VAF.races - tp_data$VAF.caller)^2))
  acc_dp <- sqrt(mean((tp_data$DP.races - tp_data$DP.caller)^2))
  acc_nv <- sqrt(mean((tp_data$NV.races - tp_data$NV.caller)^2))
  ### PLOTS ###
  p <- merged_df %>% ggplot(aes(x=VAF.races,y=VAF.caller, color=classes))+
    geom_point(na.rm = T) +
    facet_wrap(~classes) +
    my_ggplot_theme()+
    labs(x = "VAF races", y = paste0("VAF ",caller))+
    geom_miss_point() +
    ggtitle(paste0("Comparison between rRACES VAF and ",caller," VAF"),subtitle = paste0("Sensitivity: ", round(sensitivity, digits = 3), "\nPrecision: ",
                                                                                         round(precision,digits=3),
                                                                                         "\nVAF RMSE: ", round(acc_vaf,3)))

  p_dp <- merged_df %>% ggplot(aes(x=DP.races,y=DP.caller, color=classes))+
    geom_point(na.rm = T) +
    facet_wrap(~classes) +
    my_ggplot_theme()+
    labs(x = "DP races", y = paste0("DP ",caller))+
    geom_miss_point() +
    ggtitle(paste0("Comparison between rRACES DP and ",caller," DP"),subtitle = paste0("Sensitivity: ", round(sensitivity, digits = 3), "\nPrecision: ",
                                                                                        round(precision,digits=3),
                                                                                         "\nDP RMSE: ", round(acc_dp,3)))
  p_nv <- merged_df %>% ggplot(aes(x=NV.races,y=NV.caller, color=classes))+
    geom_point(na.rm = T) +
    facet_wrap(~classes) +
    my_ggplot_theme()+
    labs(x = "NV races", y = paste0("NV ",caller))+
    geom_miss_point() +
    ggtitle(paste0("Comparison between rRACES NV and ",caller," NV"),subtitle = paste0("Sensitivity: ", round(sensitivity, digits = 3), "\nPrecision: ",
                                                                                         round(precision,digits=3),
                                                                                         "\nNV RMSE: ", round(acc_nv,3)))


  filter_counts <- merged_df %>%
    dplyr::filter(type=="TP") %>%
    count(FILTER) %>%
    arrange(desc(n)) %>%
    head(top_filter)

  # Create the bar plot
  p_filter_tp <- ggplot(filter_counts, aes(x = reorder(FILTER, n), y = n)) +
    geom_bar(stat = "identity") +
    labs(x = paste0(caller," filter"), title = "Most frequent filtering label") +
    my_ggplot_theme() + coord_flip()
  pl <- plot_vaf_races+ plot_vaf_caller+p+p_filter_tp+p_dp+p_nv + plot_layout(design ='AAABBB\nCCCDDD\nCCCDDD\nEEEFFF\nEEEFFF')
  return(pl)

}


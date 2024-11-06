rm(list = ls())
library(rRACES)
library(dplyr)
#library(patchwork)
#library(ggplot2)


signatures_palette <- function(phylo_forest,seed){
  ref_path <- phylo_forest$get_reference_path()
  SBS_table_path <- gsub(pattern = "reference.fasta",replacement = "SBS_signatures.txt",x = ref_path)
  IDS_table_path <- gsub(pattern = "reference.fasta",replacement = "indel_signatures.txt",x = ref_path)
  SBS_table <-  read.csv(SBS_table_path,header=T,sep="\t")
  SBS_sign <- colnames(SBS_table)
  IDS_table <-  read.csv(IDS_table_path,header=T,sep="\t")
  IDS_sign <- colnames(IDS_table)
  sigs <- c(SBS_sign,IDS_sign)
  set.seed(seed)
  return(Polychrome::createPalette(length(sigs), c("#6B8E23","#4169E1"), M=1000,
                                   target="normal", range=c(15,80)) %>%
           setNames(sigs))
}

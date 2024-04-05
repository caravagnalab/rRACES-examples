rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils.R")

seed <- 12345
set.seed(seed)

phylo_forest <- rRACES::load_phylogenetic_forest("data/phylo_forest.sff")

# Simulate sequencing ####
seq_results <- simulate_seq(phylo_forest, coverage = 80)

seq_to_long <- function(seq_results) {
  # get names of samples
  sample_names <- strsplit(colnames(seq_results)[grepl("VAF", colnames(seq_results), fixed = TRUE)], ".VAF") %>% unlist()

  sn <- sample_names[1]
  seq_df <- lapply(sample_names, function(sn) {
    cc <- c("chromosome", "chr_pos", "ref", "alt", colnames(seq_results)[grepl(sn, colnames(seq_results), fixed = T)])

    seq_results[, cc] %>%
      `colnames<-`(c("chromosome", "chr_pos", "ref", "alt", "occurences", "coverage", "VAF")) %>%
      dplyr::mutate(sample_name = sn)
  }) %>% do.call("bind_rows", .)

  seq_df %>%
    dplyr::rename(chr = chromosome, from = chr_pos, DP = coverage, NV = occurences, ALT = alt) %>%
    dplyr::mutate(to = from) %>%
    dplyr::select(chr, from, to, ALT, NV, DP, VAF, sample_name)
}

saveRDS(seq_to_long(seq_results, "data/sequencing.rds"))

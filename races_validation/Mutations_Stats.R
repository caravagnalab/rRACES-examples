rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils/compare_snvs.R")
source("utils/vcf_parser.R")
seed <- 12345
set.seed(seed)

phylo_forest <- load_phylogenetic_forest("../SPN01/data/phylo_forest.sff")
plot <- mutation_calling_accuracy(seq_races = "/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/data/seq_results_80X_with_error_paired.rds",
				  vcf_file= "/orfeo/cephfs/scratch/cdslab/ggandolfi/races/fastq_new_version/results_SPN01_all/variant_calling/mutect2/Sample_1_vs_normal_sample/Sample_1_vs_normal_sample.mutect2.filtered.vcf.gz",
				  caller="mutect2",
				  sample_id="SPN01_Sample_1",
				  top_filter=20)
ggsave(filename = "plot_mutect2.png", plot=plot, dpi=300, width = 16, height = 16)
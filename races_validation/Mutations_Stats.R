rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils/compare_snvs.R")
source("utils/vcf_parser.R")
seed <- 12345
set.seed(seed)

#phylo_forest <- load_phylogenetic_forest("../SPN01/data/phylo_forest.sff")
seq_races = "/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/data/seq_results_100X_basic_error_350_3tumor.rds"
mutect2_dir = "/orfeo/cephfs/scratch/cdslab/ggandolfi/races/COMP_TIME/results_SPN01_MS/variant_calling/mutect2/"
patientID = "SPN01"
samples = c("Sample_1","Sample_2","Sample_3")
plot_list = list()

for (i in samples){
	vcf_file= paste0(mutect2_dir,i,"_vs_normal_sample/",
			 i,"_vs_normal_sample.mutect2.filtered.vcf.gz")
	sample_id=paste(patientID,i,sep="_")
        plot <- mutation_calling_accuracy(seq_races, vcf_file,
					  caller="mutect2",
					  sample_id,20)
	ggsave(filename = paste0("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/plots/plot_mutect2_",i,".png"), plot=plot, dpi=300, width = 16, height = 16)
}




rm(list = ls())
suppressPackageStartupMessages({
  library(ProCESS)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(optparse)
})
source("utils/compare_snvs.R")
source("utils/vcf_parser.R")
seed <- 12345
set.seed(seed)




# Define command-line options
option_list <- list(
  make_option(c("--seq_res"), type = "character", default = NULL, 
              help = "Path to the seq_results_final.rds file", metavar = "character"),
  make_option(c("--variant_calling_dir"), type = "character", default = NULL,
              help = "Path to the vcf file", metavar = "character"),
  make_option(c("--caller"), type = "character", default = NULL,
              help = "Caller used for somatic variants", metavar = "character"),
  make_option(c("--patient_id"), type = "character", default = NULL,
              help = "Patient identifier", metavar = "character"),
  make_option(c("--samples"), type = "character", default = NULL,
              help = "Comma-separated list of samples", metavar = "character"),
  make_option(c("--samples_vc"), type = "character", default = NULL,
              help = "Comma-separated list of samples for variant calling", metavar = "character")
  
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the required arguments are provided
if (is.null(opt$seq_res)) {
  stop("Error: --seq_res argument is required. Please provide the path to the seq_results_final.rds file.")
}
#if (is.null(opt$mutect2_dir)) {
#  stop("Error: --mutect2_dir argument is required. Please provide the path to the Mutect2 directory.")
#}

# Assign the command line inputs to variables
seq_races <- opt$seq_res
vc_dir <- opt$variant_calling_dir
samples <- strsplit(opt$samples, ",")[[1]]
samples_vc <- strsplit(opt$samples_vc, ",")[[1]]
caller <- opt$caller
patient_id <- opt$patient_id


get_vcf_name <- function(caller,variant_calling_dir,type,sample_id){
	if (caller=="strelka"){
		if (type=="snv"){
			vcf_name=paste0(vc_dir,"/",caller,"/",sample_id,"_vs_normal_sample/",
					sample_id,"_vs_normal_sample.strelka.somatic_snvs.vcf.gz")
		} else if (type=="indel"){
			vcf_name=paste0(vc_dir,"/",caller,"/",sample_id,"_vs_normal_sample/",
					sample_id,"_vs_normal_sample.strelka.somatic_indels.vcf.gz")
		}
	} else if (caller=="mutect2" & type=="snv"){
		vcf_name=paste0(vc_dir,"/",caller,"/",sample_id,"_vs_normal_sample/",
				sample_id,"_vs_normal_sample.mutect2.filtered.vcf.gz")
	}
        return(vcf_name)
}

#phylo_forest <- load_phylogenetic_forest("../SPN01/data/phylo_forest.sff")
#seq_races = "/fast/cdslab/ggandolfi/SPN01/tumor_purity/50X_0.3p/data/seq_results_final.rds"
#mutect2_dir = "/orfeo/cephfs/scratch/cdslab/ggandolfi/races/COMP_TIME/results_SPN01_MS/variant_calling/mutect2/"
#strelka_dir = "/orfeo/cephfs/scratch/cdslab/shared/races/sarek_SPN01/purity_ProCESS/results_SPN01_50X/variant_calling/strelka/"
#patientID = "SPN01"
#samples = c("Sample_A_03p","Sample_B_03p","Sample_C_03p")
plot_list = list()
#strelka/Sample_A_03p_vs_normal_sample/Sample_A_03p_vs_normal_sample.strelka.somatic_snvs.vcf.gz
for (i in seq_along(samples_vc)){
	vcf_file <- get_vcf_name(caller,vc_dir,"snv",samples_vc[i])
	sample_id=paste(patient_id,samples_vc[i],sep="_")
	print(sample_id)
        plot <- mutation_calling_accuracy(seq_races, vcf_file,
					  caller=caller,
					  sample_id=samples[i],20,
	                                  sample_id_vcf=sample_id)
	ggsave(filename = paste0(caller,samples_vc[i],".png"), plot=plot, dpi=300, width = 16, height = 16)
	#ggsave(filename = paste0("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/ProCESS-examples/SPN01/plots/plot_mutect2_",i,".png"), plot=plot, dpi=300, width = 16, height = 16)
}



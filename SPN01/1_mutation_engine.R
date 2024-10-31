rm(list = ls())
library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(ggrepel)
dir <- getwd()
set.seed(12345)
forest <- load_samples_forest("samples_forest_smaller.sff")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/new_simulation/sim6")
#source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/plotting/spn_blueprint/signature_palette.R")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COAD",
			   tumour_study = "US")

## Table with info ##
#
rownames <- c("Clone 1", "Clone 2", "Clone 3", "Clone 4")
colnames <- c("mu_SNV", "mu_indels", "mu_CNA")

# Create an empty dataframe with predefined row and column names
df <- data.frame(matrix(ncol = length(colnames), nrow = length(rownames)))
colnames(df) <- colnames
rownames(df) <- rownames
df <- df %>%
	mutate(mu_SNV=c(1e-8,1e-8,1e-8,1e-8),
	       mu_CNA=c(0,1e-10, 1e-12, 1e-12),
	       mu_indels=c(0,0,0,0))
tbl <- ggtexttable(df, rows = NULL, theme = ttheme("blank",base_size = 8)) %>%
          tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>%
          tab_add_hline(at.row = 4, row.side = "bottom", linewidth = 3, linetype = 1)
          #tab_add_footnote(text = "*Values referring to end of simulation", size = 8, face = "italic")

mu_SNV = 1e-8
mu_CNA = 1e-10
##112707518-112846239 
CNA_Clone2 = rRACES::CNA(type = "D", "5",
		 chr_pos = 107707518, len = 1e7,allele = 0)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "Clone 1",
		    passenger_rates = c(SNV = mu_SNV, CNA = 0),drivers = list(list("APC R1450*", allele = 1)))
m_engine$add_mutant(mutant_name = "Clone 2",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),drivers = list(CNA_Clone2))
mu_SNV = 1e-8
mu_CNA = 1e-12
m_engine$add_mutant(mutant_name = "Clone 3",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),drivers = list("TP53 R175H"))

m_engine$add_mutant(mutant_name = "Clone 4",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA),drivers = list(WGD))


# Mutational signatures
m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.40,SBS5 = 0.20,
		      SBS17a = 0.15,SBS17b = 0.15,ID1 = 0.60,ID2 = 0.40,SBS13 = 0.1))
print("Mutation engine created")
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)

#phylo_forest$save(paste0(dir,"/","phylo_forest_smaller.sff"))

setwd(dir)
########## Plotting ##########
## stick forest
labels <- get_relevant_branches(forest)
sticks <- plot_sticks(forest, labels)
## Exposure timeline

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

exp_timeline <- plot_exposure_timeline(phylo_forest)+
	ggplot2::scale_colour_manual(name = "Signatures", values = signatures_palette(phylo_forest,55)) 
## Muller plot
muller <- readRDS("muller_plot.rds")+
	theme(legend.position='right')+
	CNAqc:::my_ggplot_theme()
##
annot_forest <- plot_forest(forest) %>%
    annotate_forest(phylo_forest,
                    samples = T,
                    MRCAs = T,
                    exposures = T,
                    drivers=T,
                    add_driver_label = F)+
    ggplot2::scale_colour_manual(values = signatures_palette(phylo_forest,55))
page2 <- "AAACCC\nBBBCCC\nDDDDDD\nDDDDDD"
page2 <- "AAAAA\nBBBBB\nCCCCC\nDDDDD\nDDDDD\nEEEEE\nEEEEE"
part1 <- ggplot()+ggtitle("Muller plot")
part2 <- ggplot()+ggtitle("Exposure plot")
part3 <- ggplot()+ggtitle("Forest with sticks")
part4 <- ggplot()+ggtitle("Forest with exposure")
page2 <-wrap_plots(list(muller,exp_timeline,tbl,sticks,annot_forest),design = page2)+
	patchwork::plot_annotation(subtitle = "Exposure plots")
saveRDS(page2,"page2.rds")
ggsave("page2.pdf",plot=page2,width = 250, height = 320, units = "mm", dpi = 300)

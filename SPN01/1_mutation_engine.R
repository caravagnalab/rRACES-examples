rm(list = ls())
library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)
dir <- getwd()
set.seed(12345)
forest <- load_samples_forest("samples_forest.sff")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/new_simulation/sim6")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COAD",
			   tumour_study = "US")

mu_SNV = 1e-8
mu_CNA = 1e-10
##112707518-112846239 
CNA_Clone2 = CNA(type = "D", "5",
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

phylo_forest$save(paste0(dir,"/","phylo_forest.sff"))

setwd(dir)
########## Plotting ##########
## stick forest
labels <- get_relevant_branches(forest)
sticks <- plot_sticks(forest, labels)
## Exposure timeline
exp_timeline <- plot_exposure_timeline(phylo_forest)
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
                    add_driver_label = F)
page2 <- "AAACCC\nBBBCCC\nDDDDDD\nDDDDDD"
page2 <- "AAAAA\nBBBBB\nCCCCC\nCCCCC\nDDDDD\nDDDDD"
part1 <- ggplot()+ggtitle("Muller plot")
part2 <- ggplot()+ggtitle("Exposure plot")
part3 <- ggplot()+ggtitle("Forest with sticks")
part4 <- ggplot()+ggtitle("Forest with exposure")
page2 <-wrap_plots(list(muller,exp_timeline,sticks,annot_forest),design = page2)+
	patchwork::plot_annotation(subtitle = "Exposure plots")
saveRDS(page2,"page2.rds")
ggsave("page2.png",plot=page2,width = 250, height = 320, units = "mm", dpi = 300)

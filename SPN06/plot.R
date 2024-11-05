

rm(list = ls()) # clears objects from the workspace

#devtools::install_github("caravagnalab/rRACES")
#devtools::install_github("caravagnalab/CNAqc")

# load the required packages
#-------------------------------------------------------------------------------
library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CNAqc)
library(ggrepel)


# load the auxiliary functions
#-------------------------------------------------------------------------------
source("/u/cdslab/ahaghighi/scratch/packages/rRACES-examples/plotting/spn_blueprint/utils.R")
source("/u/cdslab/ahaghighi/scratch/packages/rRACES-examples/plotting/spn_blueprint/plot_genome_wide.R")

# load the sequence and phylogenetic forest objects
#-------------------------------------------------------------------------------
seq_results <- readRDS("data/SPN06_seq_80X.rds")
phylo_forest <- load_phylogenetic_forest("data/phylo_forest.sff")

samples <- phylo_forest$get_samples_info()[["name"]] %>% sort()
print(samples)


pdf("output/genome_wide_80X.pdf", width = 10, height = 12)





plots_gw <- list()

cov = 800
error_rate = 1e-3
tumour_type = "LUAD"
germline_sub = "default"

for (i in samples) {
  
  x <- races2cnaqc(
    seq_results=seq_results,
    phylo_forest=phylo_forest,
    sample_id=i,
    ref="GRCh38",
    purity=0.9
  )
  
  
  gw <- genome_wide_plots(x, seq_results, i)
  plots_gw[[i]] <- wrap_plots(gw[[1]], gw[[2]]) + 
    plot_layout(design="AAAAA\nAAAAA\nBBBBB\nBBBBB") + 
    plot_annotation(title = i, subtitle = paste0("Simulated coverage: ", cov,
                                                 "\nSimulated purity: ", x$purity,
                                                 "\nSequencing Errore rate: ", error_rate,
                                                 "\nTumour type: ", tumour_type,
                                                 "\nGermline subject :", germline_sub))
  
  print(plots_gw[[i]])
  seg <- wrap_plots(gw[[3]])
  print(seg)
  
}



dev.off()



#-------------------------------------------------------------------------------


# ggsave(filename = "output/genome_wide_80X.pdf", plot = p, width = 16, height = 10, dpi = 300)


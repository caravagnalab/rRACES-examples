

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
seq_results <- readRDS("data/SPN06_seq_200X.rds")
phylo_forest <- load_phylogenetic_forest("data/phylo_forest.sff")

samples <- phylo_forest$get_samples_info()[["name"]] %>% sort()
print(samples)


pdf("output/genome_wide_200X.pdf", width = 10, height = 12)

plots_gw <- list()

cov = 200
error_rate = 1e-3
tumour_type = "LUAD"
germline_sub = "default"

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
expousures_table <- phylo_forest$get_exposures() %>% dplyr::rename(causes=signature)
s_seq <- seq_results %>% filter(classes!="germinal")
s_seq_long <- s_seq %>% rRACES::seq_to_long()

for (c in unique(seq_results$chr)) {
  
  colors <- get_classes_colors(unique(s_seq$classes))
  
  plot_vaf <- s_seq_long %>%
    filter(chr==c) %>%
    filter(VAF>=0.02) %>%
    ggplot(aes(x=VAF,fill=classes)) +geom_histogram(binwidth = 0.01) +
    xlim(c(0,1))+
    facet_wrap(~ sample_name, ncol=1,scales = "free")+
    ggplot2::ggtitle(label = paste0("Chromosome ",c)) +
    CNAqc:::my_ggplot_theme()+
    ggplot2::scale_fill_manual(values =colors)
  
  p_marg <- plot_VAF_marginals(s_seq, chromosomes = c, samples = samples, labels = s_seq["classes"])
  
  if (length(p_marg)<10) {
    n <- length(p_marg)+1
    for (i in c(n:10)){
      p_marg[[i]]<-ggplot2::ggplot()
    }
  }
  
  layout <- 'AAABB\nAAABB\nAAABB\nAAABB\nAAABB'
  p <- wrap_plots(p_marg, design='AAFF\nAAFF\nBBGG\nBBGG\nCCHH\nCCHH\nDDII\nDDII\nEEJJ\nEEJJ') / plot_vaf + 
    plot_layout(guides = 'collect', design = layout) + plot_annotation(title = (paste("Chromosome", c))) & theme(legend.position = 'bottom')

  print(p)
}


#-------------------------------------------------------------------------------
exposure_prop <- phylo_forest$get_exposures() %>% 
  ggplot( aes(fill=signature, y=exposure,x=type) ) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  facet_wrap(~time)+
  CNAqc:::my_ggplot_theme()

plot_vaf <- s_seq_long %>%
  filter(VAF>=0.02) %>%
  filter(grepl("SBS|ID", causes)) %>%
  filter(!grepl("error",causes)) %>%
  inner_join(expousures_table,by="causes") %>%
  ggplot(aes(x=VAF, fill= causes,alpha=exposure)) +geom_histogram(binwidth = 0.01) +
  xlim(c(0,1))+
  facet_grid(sample_name ~ causes, scales = "free_y")+
  ggplot2::scale_alpha_continuous(range = c(0.1, 1)) + 
  ggplot2::theme_bw() +
  CNAqc:::my_ggplot_theme()

exposure_prop + plot_vaf + plot_layout(design="A\nA\n#\nB\nB\nB") + 
  plot_annotation(title = "Mutational signatures")


dev.off()


#-------------------------------------------------------------------------------


# ggsave(filename = "output/genome_wide_80X.pdf", plot = p, width = 16, height = 10, dpi = 300)


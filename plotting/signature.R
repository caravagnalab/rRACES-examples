library(rRACES)
library(ggplot2)
library(ggplot2)
library(ggalluvial)

phylo_forest <- rRACES::load_phylogenetic_forest('/orfeo/scratch/cdslab/shared/races/phylo_forest_smaller.sff')
forest <- rRACES::load_samples_forest('/orfeo/scratch/cdslab/shared/races/samples_forest.sff')

exposures <- rRACES:::get_exposure_ends(phylo_forest)
df <- exposures %>% tidyr::pivot_longer(c(time, end_time))


signatures_palette <- function(phylo_forest,seed){
  ref_path <- phylo_forest$get_reference_path()
  SBS_table_path <- '/orfeo/cephfs/scratch/cdslab/shared/races/GRCh38/SBS_signatures.txt'#gsub(pattern = "reference.fasta",replacement = "SBS_signatures.txt",x = ref_path)
  IDS_table_path <- '/orfeo/cephfs/scratch/cdslab/shared/races/GRCh38/indel_signatures.txt'#gsub(pattern = "reference.fasta",replacement = "indel_signatures.txt",x = ref_path)
  SBS_table <-  read.csv(SBS_table_path,header=T,sep="\t")
  SBS_sign <- colnames(SBS_table)
  IDS_table <-  read.csv(IDS_table_path,header=T,sep="\t")
  IDS_sign <- colnames(IDS_table)
  sigs <- c(SBS_sign,IDS_sign)
  set.seed(seed)
  return(Polychrome::createPalette(length(sigs), c("#6B8E23","#4169E1"), M=1000, target="normal", range=c(15,80)) %>% setNames(sigs))
}

sanky <- ggplot(df, aes(x = factor(round(value,2)), stratum = signature, alluvium = signature, y = exposure)) +
  geom_alluvium(aes(fill = signature), width = 0.1) +
  scale_fill_manual('Signature', values = signatures_palette(phylo_forest,55))   +
  geom_stratum(width = 0.1, alpha = .25) +
  scale_x_discrete(expand = c(.05, .05)) +
  labs(
    x = "Generation",
    y = "Exposure") +
  theme_minimal() +
  facet_wrap(.~type, nrow = 2) 

muller_sanky <- sanky + (muller+theme(legend.position = 'right')) + plot_layout(nrow = 2, design = 'A\nA\nB') 
saveRDS(object = muller_sanky, file = '/orfeo/LTS/LADE/LT_storage/lvaleriani/races/sanky.RDS')

toplot <- lapply(min(df$value):max(df$value), function(t){
  df %>% 
    filter(value <= t) %>% 
    filter(value == max(value)) %>% 
    mutate(value=t)
}) %>% do.call(rbind, .)

tree_sign <- toplot %>% 
  group_by(value, type) %>% 
  reframe(exposure = exposure/sum(exposure), across(everything())) %>% 
  ggplot()+
  geom_bar(aes(x = value, y = exposure, fill = signature), stat = "identity", width = 1) +
  ylab('Exposure') + 
  xlab('') + 
  scale_x_reverse() +
  scale_fill_manual('Signature', values = signatures_palette(phylo_forest,55))+ 
  coord_flip() + 
  theme_minimal() + 
  facet_grid(.~type)

final_forest <- plt_forest + tree_sign
saveRDS(object = tree_sign, file = '/orfeo/LTS/LADE/LT_storage/lvaleriani/races/tree_sign.RDS')

signature_plt <- wrap_plots(muller_sanky, final_forest, nrow = 2)
ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/signature.png', signature_plt, width = 210, height = 297, units = "mm")


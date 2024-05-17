library(rRACES)
library(dplyr)
library(ggplot2)

setwd('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/')
forest <- load_samples_forest("/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/samples_forest.sff")
m_engine <- build_mutation_engine(setup_code = "GRCh38", context_sampling = 20)
plot_forest <- T

# Clone 1 
# NOTCH1p2514*fs*4 
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = 3e-8,
                                        CNA = 1e-10),
                    driver_SNVs = list("NOTCH1 FY357Y")
)

# Clone 2 
# KRAS G12D
m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV = 3e-8, 
                                        CNA = 1e-10),
                    driver_SNVs = list("KRAS G12D")
)

# Clone 3
# Unknown
m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = 3e-8, 
                                        CNA = 1e-10),
                    drivers = list(SNV("1", 59550210, "A"))
)

# Signatures
# SBS1, SBS5, IDSXX
m_engine$add_exposure(c(SBS5 = 0.5, SBS1 = 0.5))
m_engine

phylo_forest <- m_engine$place_mutations(forest, 1000)

all_SNV <- phylo_forest$get_sampled_cell_mutations() %>% as_tibble()
all_SNV %>%
  group_by(cause) %>%
  summarise(nPos = n_distinct(chr_pos))

all_SNV %>%
  group_by(class) %>%
  summarise(nPos = n_distinct(chr_pos))

all_CNV <- phylo_forest$get_sampled_cell_CNAs() %>% as_tibble()
all_CNV %>%
  group_by(class) %>%
  summarise(nPos = n_distinct(begin))

phylo_forest$save("/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/results/phylo_forest.sff")

if (plot_forest == T){
  annot_forest <- plot_forest(forest) %>%
    annotate_forest(phylo_forest, exposures = F)
  ggsave('/orfeo/LTS/LADE/LT_storage/lvaleriani/races/SPN03/plots/SPN03_ann_forest.png', plot = annot_forest,   height = 12, width = 10, dpi = 300, units = 'in')
}


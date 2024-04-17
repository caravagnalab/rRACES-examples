library(rRACES)
library(dplyr)
library(ggplot2)

forest <- load_samples_forest("./results/samples_forest.sff")
m_engine <- build_mutation_engine(setup_code = "GRCh38", context_sampling = 20)
plot_forest <- F

# Clone 1 
# NOTCH1p2514*fs*4 
# 2 bp frameshift deletion
# GRCh37.p13 chr:9	NC_000009.11:g.139390649_139390650del
# GRCh38.p14 chr:9	NC_000009.12:g.136496197_136496198del
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = 3e-8,
                                        CNA = 1e-10),
                    driver_SNVs = c(SNV("9", 136496197, "C"))
)

# Clone 2 
# KRAS G12D
# GRCh37 chr:12 start:25398284 end:25398284 C>T
m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV = 3e-8, 
                                        CNA = 1e-10),
                    driver_SNVs = c(SNV("12", 25398284, "T"))
)

# Clone 3
# Unknown
m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = 3e-8, 
                                        CNA = 1e-10),
                    driver_SNVs = c(SNV("1", 59550210, "A"))
)

# Signatures
# SBS1, SBS5, IDSXX
m_engine$add_exposure(c(SBS5 = 0.5, SBS1 = 0.5))
m_engine

phylo_forest <- m_engine$place_mutations(forest, 1000)

all_SNV <- phylo_forest$get_sampled_cell_SNVs() %>% as_tibble()
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

phylo_forest$save("./results/phylo_forest.sff")

if (plot_forest == T){
  annot_forest <- plot_forest(forest) %>%
    annotate_forest(phylo_forest, exposures = F)
  ggsave('./plots/SPN03_ann_forest.png', plot = annot_forest,   height = 12, width = 10, dpi = 300, units = 'in')
}


library(rRACES)
library(tidyverse)

phylo_forest <- load_phylogenetic_forest("./results/phylo_forest.sff")
all_SNV <- phylo_forest$get_sampled_cell_SNVs() %>% as_tibble()

all_SNV %>%
  group_by(cause) %>%
  summarise(nPos = n_distinct(chr_pos))

all_SNV %>%
  group_by(class) %>%
  summarise(nPos = n_distinct(chr_pos))

cov <- 80
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

seq_results <- parallel::mclapply(chromosomes, function(c) {
        simulate_seq(phylo_forest, coverage = cov, chromosomes = c, write_SAM = FALSE)
}, mc.cores = parallel::detectCores()) %>% do.call("bind_rows", .)

saveRDS(seq_results, file = paste0('./results/seq_', cov, 'X.RDS'))

# drivers validation
library(tidyverse)
library(ProCESS)

# getters
source('getters/process_getters.R')
source('getters/sarek_getters.R')

samples = get_sample_names(spn = 'SPN03')
phylo_forest = get_phylo_forest('SPN03', base_path = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT')
phylo_forest = load_phylogenetic_forest(phylo_forest)

drivers = get_drivers_results(spn = 'SPN03', samples = samples, purity = '0.9', coverage = '50x', callers = 'mutect2_ascat', cohort = 'SCOUT')
drivers_tumourevo = lapply(drivers, function(x) {
  x[[1]]$mutations %>% 
    dplyr::filter(is_driver)
})

# process_drivers = phylo_forest$get_driver_mutations()
process_seq_res = get_process_seq(spn = 'SPN03', coverage = '50x', purity = 0.9)

process_drivers = process_seq_res %>% 
  dplyr::filter(classes == 'driver')

drivers_tumourevo = lapply(names(drivers_tumourevo), function(x) {
  drivers_tumourevo[[x]] %>% 
    dplyr::mutate(sample = x) %>% 
    dplyr::select(chr, from, to, ref, alt, NV, DP, VAF, driver_label, sample)
})
names(drivers_tumourevo) = lapply(drivers_tumourevo, function(x) {
  x$sample %>% unique
}) %>% unlist

driver_comparison = lapply(names(drivers_tumourevo), function(x) {
  
  te = drivers_tumourevo[[x]]
  te = te %>% 
    dplyr::rename(tumourevo_NV = NV) %>% 
    dplyr::rename(tumourevo_DP = DP) %>% 
    dplyr::rename(tumourevo_VAF = VAF)
  
  pr = process_drivers %>% 
    dplyr::select(chr, chr_pos, ref, alt, causes, starts_with(x)) %>% 
    dplyr::mutate(chr = paste0('chr', chr))
  
  colnames(pr) = c('classes', 'chr', 'from', 'ref', 'alt', 'causes', 'process_NV', 'process_DP', 'process_VAF')
  
  full_join(te, pr)
  
})

names(driver_comparison) = names(drivers_tumourevo)




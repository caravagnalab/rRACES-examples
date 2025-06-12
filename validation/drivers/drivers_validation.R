# validation of driver calls
get_drivers_results = function(path = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT', 
                               spn, 
                               samples, 
                               purity, 
                               coverage, 
                               callers, 
                               cohort) {
  
  drivers_path = paste0(path, '/', spn, '/tumourevo/', coverage, '_', purity, 'p_', callers, '/driver_annotation/annotate_driver/', cohort, '/', spn)
  drivers_path = lapply(samples, function(x) {
    list.files(drivers_path, pattern = x, full.names = T, recursive = T)
  })
  
  drivers = lapply(drivers_path, readRDS)
  names(drivers) = samples
  return(drivers)
}

# get process sequencing
get_process_seq = function(path = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT', 
                           spn, 
                           coverage, 
                           purity
) {
  seq_path = paste0(path, '/', spn, '/sequencing/tumour/purity_', purity, '/data/mutations/seq_results_muts_merged_coverage_', coverage, '.rds')
  seq_obj = readRDS(seq_path)
  return(seq_obj)
}

# plot comparison 

plot_comparison_drivers = function(x) {
  x = x %>% 
    dplyr::mutate(causes = ifelse(is.na(causes), 'Tumourevo annotation', causes)) %>% 
    dplyr::select(-c(to)) %>% 
    dplyr::mutate(across(matches('NV|DP|VAF'), ~ ifelse(is.na(.), 0, .)))
  
  
  x %>% 
    ggplot(aes(x = tumourevo_VAF, y = process_VAF, colour =  classes, label = driver_label)) + 
    geom_point() + 
    facet_wrap(vars(causes), scales = 'free') + 
    theme_bw() + 
    xlim(c(-0.01, 1.01)) + 
    ylim(c(-0.01, 1.01)) + 
    ggrepel::geom_text_repel(max.overlaps = 1)
  
}



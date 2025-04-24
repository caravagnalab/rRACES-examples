get_phylo_forest = function(sample, local = F){
  if (local){philo_dir = '~/dati_Orfeo/shared/SCOUT/'}
  phylo_forest_path = paste0('~/dati_Orfeo/shared/SCOUT/',sample,'/races')
  rRACES::load_phylogenetic_forest(phylo_forest_path)
}
phylo_forest_path = get_phylo_forest('SPN01')
sample_id <- opt$sample_id
segment_file<-opt$segment_file
overlap_threshold<-opt$overlap_threshold
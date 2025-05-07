library(ProCESS)
library(dplyr)
library(optparse)

option_list <- list( 
  make_option(c("-i", "--input"), type="character", default='/orfeo/cephfs/scratch/cdslab/shared/SCOUT', help="path to input data"),
  make_option(c("-s", "--SPN"), type="character", default='SPN03', help="SPN name"),
  make_option(c("-c", "--coverage"), type="character", default='50', help="coverage"),
  make_option(c("-p", "--purity"), type="character", default='0.6', help="purity")
  )

param <- parse_args(OptionParser(option_list=option_list))
dir <- param$input
spn <- param$SPN
coverage <- param$coverage
purity <- param$purity
out <- paste0(dir, '/', spn, '/process/cna_data') 


seq_to_long <- function(seq_results) {
  sample_names <- strsplit(colnames(seq_results)[grepl(".VAF", colnames(seq_results), fixed = TRUE)], ".VAF") %>% unlist()
  
  seq_df <- lapply(sample_names, function(sn) {
    cc <- c("chr", "chr_pos", "ref", "alt", "causes", "classes", colnames(seq_results)[grepl(paste0(sn, "."), colnames(seq_results), fixed = TRUE)])
    seq_results[, cc] %>%
      `colnames<-`(c("chr", "chr_pos", "ref", "alt", "causes", "classes", "occurences", "coverage", "VAF")) %>%
      dplyr::mutate(sample_name = sn)
  }) %>% do.call("bind_rows", .)
  
  seq_df %>%
    dplyr::rename(chr = chr, from = chr_pos, DP = coverage, NV = occurences, ALT = alt) %>%
    dplyr::mutate(to = from)
}

print('Read rds')
# read rds
rds <- readRDS(paste0(dir, '/', spn, '/sequencing/tumour/purity_', purity, '/data/mutations/seq_results_muts_merged_coverage_', coverage, 'x.rds'))
print('Done')

print('Transform to long format')
# transform to long
long_rds <- seq_to_long(rds)
long_rds_sample <- lapply(unique(long_rds$sample_name), FUN = function(sample){
  long_rds %>%
    filter(sample_name == sample)
})
names(long_rds_sample) <- unique(long_rds$sample_name)
print('Done')

print('Read cnas data')
# read cnas
files_cna <- list.files(paste0(dir, '/', spn, '/process/cna_data'), full.names = T, pattern = 'cna.rds')
sample_names <- sapply(files_cna, function(path) {
  base_name <- basename(path)
  sub("_cna.rds$", "", base_name)
})
cna_seg <- lapply(files_cna, FUN = function(f){
  base_name <- basename(f)
  sample <- sub("_cna.rds$", "", base_name)
  readRDS(f) %>%
    dplyr::mutate(CN_type = ifelse(ratio < 0.9 & ratio > 0.1, 'sub-clonal', 'clonal'),
           CN = paste(major, minor, sep = ':'),
          seg_id = paste(chr,begin,end, sep = ':'),
          sample = sample)
})
names(cna_seg) <- sample_names
print('Done')

print('Map mutations to cnas')
muts_cn <- list()
for (sample in names(cna_seg)){
    print(sample)
    tmp_cna <- cna_seg[[sample]]

    muts_cn[[sample]] <- long_rds_sample[[sample]] %>%
      filter(classes != "germinal") %>% 
      inner_join(tmp_cna, by = c("chr"), relationship = "many-to-many") %>%
      filter(from >= begin & to <= end)
}
print('Done')

print('Save rds')
saveRDS(object = muts_cn, file = paste0(out, '/cna_muts_purity_',param$pur, '_coverage_', param$cov,'x.rds'))
print('Done')


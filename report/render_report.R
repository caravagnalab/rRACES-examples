library(dplyr)
library(optparse)

option_list <- list( 
  make_option(c("-w", "--workdir"), type="character", default='/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/report', help="path to ProCESS-examples/report directory"),
  make_option(c("-i", "--input"), type="character", default='/orfeo/cephfs/scratch/cdslab/shared/SCOUT', help="path to input directory"),
  make_option(c("-s", "--SPN"), type="character", default='SPN03', help="SPN name"),
  make_option(c("-o", "--output"), type="character", default='/orfeo/cephfs/scratch/cdslab/shared/SCOUT', help="path to output directory")
)

COV = c(50, 100, 150, 200)
PUR = c(0.3, 0.6, 0.9)

opt <- parse_args(OptionParser(option_list=option_list))
spn <- opt$SPN

dir.create(paste0(opt$output,'/', spn, '/report/'), recursive = T, showWarnings = F)
setwd(opt$workdir)

sim <- paste0(opt$input, '/', spn, '/process')
sample_forest <- paste0(opt$input, '/', spn, '/process/sample_forest.sff')
phylo_forest <- paste0(opt$input, '/', spn, '/process/phylo_forest.sff')
cna_dir <- paste0(opt$input, '/', spn, '/process/cna_data')
seq_res <- paste0(opt$input, '/', spn, '/sequencing')


print('Rendering SPN report')
params <- list(
  spn = spn,
  files = list(
    base_dir = sim,
    sample_forest = sample_forest,
    phylo_forest = phylo_forest,
    seq_res = seq_res
  ),
  cna_dir = cna_dir
)
rmarkdown::render("Report.Rmd", 
                  params = params, 
                  output_file = paste0(opt$output, '/', spn, '/report/', spn, '.html'))

print('Done')

# for (cov in COV){
#   for (pur in PUR){
#     file_name <- paste0('purity_', pur,'_coverage_', cov, 'x_', spn, '.yml')
#     print(paste0('Rendering SPN sequencing report: purity ', pur, ', coverage ', cov))
#     
#     config <- list(
#       spn = spn,
#       sequencing = list(
#         coverage = as.integer(cov),
#         purity = pur
#       ),
#       files = list(
#         sim = sim,
#         sample_forest = sample_forest,
#         phylo_forest = phylo_forest,
#         seq_res = seq_res
#       ),
#       cna_dir = cna_dir
#     )
#     rmarkdown::render("Report_Seq.Rmd", 
#                       params = config, 
#                       output_file = paste0(param$output,'/', spn, '/report/', spn, '_purity_', pur, '_coverage_', cov, 'x.html'))
#     print('Done')
#   }
# }
# 


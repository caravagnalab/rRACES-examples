library(dplyr)

COV = c(50, 100, 150, 200)
PUR = c(0.3, 0.6, 0.9)

########################
# change with SPN name
spn <- "SPN03"
########################
workdir  <- "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/report"
output <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
input <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"

dir.create(paste0(output,'/', spn, '/report/'), recursive = T, showWarnings = F)
setwd(workdir)

print('Rendering SPN report:')
params <- list(
  spn = spn,
  files = list(
    base_dir = paste0(input, '/', spn, '/process'),
    sample_forest = paste0(input, '/', spn, '/process/sample_forest.sff'),
    phylo_forest = paste0(input, '/', spn, '/process/phylo_forest.sff'),
    seq_res = paste0(input, '/', spn, '/sequencing')
  ),
  cna_dir = paste0(input, '/', spn, '/process/cna_data')
)
rmarkdown::render("Report.Rmd", 
                  params = params, 
                  output_file = paste0(output, '/', spn, '/report/', spn, '.html'))

print('Done')

for (cov in COV){
  for (pur in PUR){
    file_name <- paste0('purity_', pur,'_coverage_', cov, 'x_', spn, '.yml')
    print(paste0('Rendering SPN sequencing report: purity ', pur, ', coverage ', cov))

    params <- list(
      spn = paste0(input, '/', spn, '/process'),
      sequencing = list(
        coverage = as.integer(cov),
        purity = pur
      ),
      files = list(
        sim =  paste0(input, '/', spn, '/process'),
        sample_forest = paste0(input, '/', spn, '/process/sample_forest.sff'),
        phylo_forest = paste0(input, '/', spn, '/process/phylo_forest.sff'),
        seq_res = paste0(input, '/', spn, '/sequencing')
      ),
      cna_dir = paste0(input, '/', spn, '/process/cna_data')
    )
    rmarkdown::render("Report_Seq.Rmd",
                          params = params,
                      output_file = paste0(output,'/', spn, '/report/', spn, '_purity_', pur, '_coverage_', cov, 'x.html'))
    print('Done!')
  }
}

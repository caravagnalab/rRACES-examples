library(dplyr)

# TO CHANGE #
########################
# change with SPN name
spn <- "SPN01"
# change with current workdir
workdir  <- "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/"
########################

output <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
input <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"

COV = c(50)
PUR = c(0.6)
#NORMAL_COV = 30

dir.create(paste0(output,'/', spn, '/validation/report/'), recursive = T, showWarnings = F)
setwd(workdir)

for (cov in COV){
  for (pur in PUR){
    comb <- paste0(cov,"x_",pur,"p")
    params <- list(
      spn = spn,
      comb = comb
    )
    print(paste0('Rendering SPN validation report: purity ', pur, ', coverage ', cov))
    rmarkdown::render("Validation_Report.Rmd",
                      params = params,
                      output_file = paste0(output,'/', spn, '/validation/report/', spn, '_purity_', pur, '_coverage_', cov, 'x.html'))
  }
}

print('Done')

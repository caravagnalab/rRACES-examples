coverages <- c(50,100)
purities <- c(0.3,0.9)
spn_id <- "SPN01"
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/REPO_UPDATED/rRACES-examples/report/")
for (coverage in coverages) {
  for (purity in purities){
    params_yaml <- paste0("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/REPO_UPDATED/rRACES-examples/SCOUT/",spn_id,"/",
                          "params_",coverage,"x_",purity,"p.yml")
    params_list <- yaml::read_yaml(params_yaml)
    outfile <- paste0("Report_Seq_",spn_id,"_",coverage,"x_",purity,"p.html")
    rmarkdown::render(input = "Report_Seq_v1.Rmd",params = params_list,
                      output_dir = "/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/REPO_UPDATED/rRACES-examples/report/",
                      output_file = outfile) 
  }
}

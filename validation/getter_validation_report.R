get_mutations_metrics <- function(spn,comb,samples,type){
  scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
  if (type=="somatic"){
    validation_dir <- paste0(scout_dir,spn,"/validation/somatic/",spn)
    tools <- c("mutect2","strelka","freebayes")
  } else if (type=="germline"){
    validation_dir <- paste0(scout_dir,spn,"/validation/germline/report")
    tools <- c("haplotypecaller","strelka","freebayes")
  }
  
  
  all_metrics <- list()
  for (sample in samples){
    all_metrics_tool <- list()
    for (tool in tools){
      if (type=="somatic"){
        purity <- strsplit(comb,"_")[[1]][2] %>% stringr::str_replace(pattern = "p",replacement = "")
        coverage <- strsplit(comb,"_")[[1]][1] %>% stringr::str_replace(pattern = "x",replacement = "")
        somatic_indel_rds <- readRDS(paste0(validation_dir,"/",comb,"/",tool,"/",sample,"/INDEL/metrics.rds"))[["report_metrics"]] %>% 
          dplyr::mutate("Tool"=tool) %>% 
          dplyr::mutate("Purity"=purity) %>% 
          dplyr::mutate("Coverage"=coverage) %>% 
          dplyr::mutate("Sample"=sample) %>% 
          dplyr::mutate("Type"="INDEL")
        somatic_snvs_rds <- readRDS(paste0(validation_dir,"/",comb,"/",tool,"/",sample,"/SNV/metrics.rds"))[["report_metrics"]] %>% 
          dplyr::mutate("Tool"=tool) %>% 
          dplyr::mutate("Purity"=purity) %>% 
          dplyr::mutate("Coverage"=coverage) %>% 
          dplyr::mutate("Sample"=sample) %>% 
          dplyr::mutate("Type"="SNV")
        all_metrics_tool[[tool]] <- rbind(somatic_indel_rds,somatic_snvs_rds) %>% 
          dplyr::mutate(value=round(as.numeric(value),2))
      } else if (type=="germline"){
          germline_rds <- readRDS(paste0(validation_dir,"/",tool,"_normal_metrics.rds"))[["report_metrics"]] %>% 
            dplyr::mutate("Tool"=tool)
          all_metrics_tool[[tool]] <- germline_rds %>% 
            dplyr::mutate(value=round(as.numeric(value),2))
      }
    }
    all_metrics[[sample]] <- do.call("rbind",all_metrics_tool)
  }

  mutation_metrics <- do.call("rbind",all_metrics)
  return(mutation_metrics)
}


get_cna_metrics <- function(spn,comb,samples){
  scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
  validation_dir_cna <- paste0(scout_dir,spn,"/validation/cna/",spn)
  tools_cna <- c("ascat")
  all_metrics <- list()
  for (sample in samples){
    all_metrics_tool <- list()
    for (tool in tools_cna){
      purity <- strsplit(comb,"_")[[1]][2] %>% stringr::str_replace(pattern = "p",replacement = "")
      coverage <- strsplit(comb,"_")[[1]][1] %>% stringr::str_replace(pattern = "x",replacement = "")
      cna_data <- readRDS(paste0(validation_dir_cna,"/",comb,"/",tool,"/",sample,"/metrics.rds"))
      cna_metrics <- stack(cna_data[5:length(cna_data)])
      colnames(cna_metrics) <- c("value","name")
      cna_metrics <- cna_metrics %>% 
        dplyr::mutate("Tool"=tool) %>% 
        dplyr::mutate("Purity"=purity) %>% 
        dplyr::mutate("Coverage"=coverage) %>% 
        dplyr::mutate("Sample"=sample)
      all_metrics_tool[[tool]] <- cna_metrics
    }
    all_metrics[[sample]] <- do.call("rbind",all_metrics_tool)
  }
  
  cna_metrics <- do.call("rbind",all_metrics)
  return(cna_metrics)
}

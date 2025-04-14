library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(paste("Syntax error: ProCESS_merge_rds.R",
             "<num_of_lots> <SPN> <input_dir> <purity> <type>"),
       call. = FALSE)
}



lot_end <- as.double(args[1])
spn <- args[2]
input_dir <- args[3]
purity <- args[4]
type <- args[5]

if (type=="tumour"){
  muts_dir <- paste0(input_dir,"tumour/purity_",purity,"/data/mutations/")
  print(muts_dir)
  max_coverage <- 200 ## this is hard-coded now
  num_of_lots <- 40 ## this is hard-coded now
  coverage<-(max_coverage*lot_end)/num_of_lots
  data <- list()
  
  if (lot_end==10){
    rds_files <- list.files(path = muts_dir, pattern = paste0("seq_results_muts_",spn,"_"),full.names = T)[1:lot_end]
    data <- lapply(rds_files,function(x){
      readRDS(x)  %>%
        dplyr::select(-ends_with(".VAF"))
    })
  } else{
    lot_start<-lot_end-10+1
    previous_lot <- lot_end-10
    previous_coverage <- (max_coverage*previous_lot)/num_of_lots
    previous_lots <- list.files(path = muts_dir, pattern = paste0("seq_results_muts_merged_coverage_",previous_coverage),full.names = T)
    rds_files <- list.files(path = muts_dir, pattern = paste0("seq_results_muts_",spn,"_"),full.names = T)[lot_start:lot_end]
    
    rds_files_all <- c(rds_files,previous_lots)
    data <- lapply(rds_files_all,function(x){
      readRDS(x)  %>%
        dplyr::select(-ends_with(".VAF"))
    })
  }
  
  ids <- grep(pattern = "coverage",x = colnames(data[[1]]),value = T) %>% strsplit("\\.")
  print("Combining dataframes ...")
  combined_df <- bind_rows(data)
  sample_names <- sapply(ids, function(x) {
    parts <- unlist(strsplit(x, "\\."))
    paste(parts[1], parts[2], sep = ".")
  })
  
  print("Summing up NV and DP...")
  
  columns  <- colnames(combined_df)[c(7:ncol(combined_df))]
  result <- combined_df %>%
    group_by(chr, chr_pos,ref,alt,classes,causes) %>% 
    summarize(across(all_of(columns), sum, .names = "{.col}"))
  
  print("Recalculate VAF...")
  for (s in sample_names){
    col_name_DP <- paste0(s,".coverage")
    col_name_NV <- paste0(s,".occurrences")
    col_name_VAF <- paste0(s,".VAF")
    result <- result %>% 
      mutate(!!col_name_VAF := .data[[col_name_NV]] / .data[[col_name_DP]])
    print(s)
  }
  print("Saving merged rds...")
  saveRDS(result, file = paste0(muts_dir,"seq_results_muts_merged_coverage_",coverage,"x", ".rds"))
  print("Done merging!")
} else if (type=="normal"){
  muts_dir <- paste0(input_dir,"normal/purity_1/data/mutations/")
  rds_files <- list.files(path = muts_dir, pattern = paste0("seq_results_muts_",spn,"_"),full.names = T)
  data <- lapply(rds_files,function(x){
    readRDS(x)  %>%
      dplyr::select(-ends_with(".VAF"))
  })
  print("Combining dataframes ...")
  combined_df <- bind_rows(data)
  
}

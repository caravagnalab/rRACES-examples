

# return corresponding matched directory based on spn, coverage, purity, vcf_caller and cna_caller
dir_getter <- function(
    spn, 
    coverage, 
    purity, 
    vcf_caller, 
    cna_caller, 
    base_path="/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
) {
  
  vcf_caller_list <- c("mutect2", "strelka")
  cna_caller_list <- c("ascat")
  
  if (!is.numeric(coverage) || coverage %% 1 != 0) {
    stop("Error: 'coverage' must be integer") # coverage : integer
  }
  if (!is.numeric(purity) || (purity > 1 || purity < 0)) {
    stop("Error: 'purity' must be a number between 0 and 1") # purity : float [0, 1]
  }
  if (!(vcf_caller %in% vcf_caller_list)) {
    stop("Error: 'vcf_caller' is not valid")
  }
  if (!(cna_caller %in% cna_caller_list)) {
    stop("Error: 'typecna_caller' is not valid")
  }
  
  # EXTEND PATH
  MAIN_PATH <- file.path(base_path, spn, "tumourevo")
  
  ptr <- paste(paste0(coverage, "x"), paste0(purity, "p"), vcf_caller, cna_caller, sep = "_")
  all_dirs <- list.dirs(path = MAIN_PATH, full.names = T, recursive = FALSE)
  matched_dir <- all_dirs[grepl(ptr, basename(all_dirs))]
  
  return(matched_dir)
}



# given a directory path as input, returns the named list of files residing in that directory
get_named_file_list <- function(dir_path) {
  files <- list.files(dir_path, full.names = TRUE)
  file_names <- basename(files)
  names(files) <- tools::file_path_sans_ext(file_names)
  return(as.list(files))
}



mobster_named_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    file_path <- input_list[[i]]
    file_name <- basename(file_path)
    
    # Extract extension
    file_ext <- tools::file_ext(file_name)
    
    # Get base name without extension
    name_without_ext <- tools::file_path_sans_ext(file_name)
    
    # Strip everything before "mobsterh" or "REPORT"
    simplified_key <- sub(".*?(mobsterh.*|REPORT.*)", "\\1", name_without_ext)
    
    # Make key unique if already present
    #if (simplified_key %in% names(output_list)) {
    simplified_key <- paste0(simplified_key, "_", file_ext)
    #}
    
    output_list[[simplified_key]] <- file_path
  }
  
  return(output_list)
}


ctree_named_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    file_path <- input_list[[i]]
    file_name <- basename(file_path)
    
    # Extract extension
    file_ext <- tools::file_ext(file_name)
    
    # Get base name without extension
    name_without_ext <- tools::file_path_sans_ext(file_name)
    
    # Strip everything before "mobsterh" or "REPORT"
    simplified_key <- sub(".*?(ctree*|REPORT.*)", "\\1", name_without_ext)
    
    # Make key unique if already present
    #if (simplified_key %in% names(output_list)) {
    simplified_key <- paste0(simplified_key, "_", file_ext)
    #}
    
    output_list[[simplified_key]] <- file_path
  }
  
  return(output_list)
}



pyclone_input_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    file_path <- input_list[[i]]
    file_name <- basename(file_path)
    
    # Remove "SCOUT_SPNxx_" prefix (xx can be any number)
    name_without_prefix <- sub("^SCOUT_SPN\\d+_", "", file_name)
    
    # Replace dots in extensions with underscores
    name_cleaned <- gsub("\\.", "_", name_without_prefix)
    
    output_list[[name_cleaned]] <- file_path
  }
  
  return(output_list)
}


viber_input_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    file_path <- input_list[[i]]
    file_name <- basename(file_path)
    
    # Remove "SCOUT_SPNxx_" prefix (xx can be any number)
    name_without_prefix <- sub("^SCOUT_SPN\\d+_", "", file_name)
    
    # Replace dots in extensions with underscores
    name_cleaned <- gsub("\\.", "_", name_without_prefix)
    
    output_list[[name_cleaned]] <- file_path
  }
  
  return(output_list)
}




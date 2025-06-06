#!/usr/bin/env Rscript

# return list of files from sarek for a given caller/purity/coverage combination
get_sarek_variant_called_files <- function(spn,
                                           sampleID,
                                           coverage,
                                           purity,
                                           variant_caller,
                                           type,
                                           basedir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/",
                                           normalID="normal_sample",
                                           patientID=NULL) {
  # testing
  if (!is.null(type)){
    if (type %in% c("tumour", "normal")) {
      stop("Options for type are tumour or normal")
    }
  }
  # input checking
  accepted_callers <- c("mutect2", "strelka", "ascat", "freebayes", "haplotypecaller", "cnvkit")
  if (!(variant_caller %in% accepted_callers)) {
    stop("Variant caller not supported. Or check spelling of variant_caller arg.")
  }
  
  # get files
  if (type == "tumour" || is.null(type)) {
    if (variant_caller %in% c("strelka", "ascat", "cnvkit", "freebayes", "haplotypecaller")) {
      sample_naming <- paste0(sampleID, "_vs_", normalID)
      path_to_files <- file.path(basedir, spn, "sarek",
                                 paste0(as.character(coverage),
                                       "x_",
                                       as.character(purity),
                                       "p"),
                                 "variant_calling",
                                 variant_caller,
                                 sample_naming)
      output_files <- list.files(path_to_files, full.names = TRUE)
    } else if (variant_caller == "mutect2") {
      sample_naming <- spn
      path_to_files <- file.path(basedir, spn, "sarek",
                                  paste0(as.character(coverage),
                                        "x_",
                                        as.character(purity),
                                        "p"),
                                  "variant_calling",
                                  variant_caller,
                                  sample_naming
                                  )
      output_files <- list.files(path_to_files, full.names = TRUE)
    } else if (variant_caller == "haplotypecaller") {
      stop("Error: No tumour files associated with caller: haplotypecaller")
    }
  } else if (type=="normal") {
    if (variant_caller %in% c("strelka", "freebayes", "haplotypecaller")) {
      sample_naming <- normalID
      path_to_files <- file.path(basedir, spn, "sarek",
                                 paste0(as.character(coverage),
                                        "x_",
                                        as.character(purity),
                                        "p"),
                                 "variant_calling",
                                 variant_caller,
                                 sample_naming)
      output_files <- list.files(path_to_files, full.names = TRUE)
    } else if (variant_caller == "mutect2") {
      stop("Error: No normal files associated with caller: mutect2")
    }
  }
  return(output_files)
}


# take a list of filenames from strelka, mutect2, or ascat and structure them 
# into a named list
parse_sarek_variant_called_files <- function(list_of_output_files) {
  # check which caller we are looking at
  if (grepl("mutect", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- vector("list")
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "filtered.vcf.gz")) {
        named_files[["vcf"]] <- list_of_output_files[i]
      } else if (endsWith(list_of_output_files[i], "filtered.vcf.gz.tbi")) {
        named_files[["tbi"]] <- list_of_output_files[i]
      }
    }
  } else if (grepl("strelka", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- vector("list")
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "snvs.vcf.gz")) {
        named_files[["snvs_vcf"]] <- list_of_output_files[i]
      } else if (endsWith(list_of_output_files[i], "snvs.vcf.gz.tbi")) {
        named_files[["snvs_tbi"]] <- list_of_output_files[i]
      } else if (endsWith(list_of_output_files[i], "indels.vcf.gz")) {
        named_files[["indels_vcf"]] <- list_of_output_files[i]
      } else if (endsWith(list_of_output_files[i], "indels.vcf.gz.tbi")) {
        named_files[["indels_tbi"]] <- list_of_output_files[i] 
      } else if (endsWith(list_of_output_files[i], "variants.vcf.gz")) {
        named_files[["variants_vcf"]] <- list_of_output_files[i] 
      } else if (endsWith(list_of_output_files[i], "variants.vcf.gz.tbi")) {
        named_files[["variants_tbi"]] <- list_of_output_files[i] 
      }
    }
  } else if (grepl("freebayes", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- vector("list")
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "freebayes.vcf.gz")) {
        named_files[["vcf"]] <- list_of_output_files[i]
      } else if (endsWith(list_of_output_files[i], "freebayes.vcf.gz.tbi")) {
        named_files[["tbi"]] <- list_of_output_files[i]
      }
    }
  } else if (grepl("haplotypecaller", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- vector("list")
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "haplotypecaller.filtered.vcf.gz")) {
        named_files[["vcf"]] <- list_of_output_files[i]
      } else if (endsWith(list_of_output_files[i], "haplotypecaller.filtered.vcf.gz.tbi")) {
        named_files[["tbi"]] <- list_of_output_files[i] 
      }
    }
  } else if (grepl("ascat", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- vector("list")
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "purityploidy.txt")) {
        named_files[["purity_ploidy"]] <- list_of_output_files[i] 
      } else if (endsWith(list_of_output_files[i], "segments.txt")) {
        named_files[["segments"]] <- list_of_output_files[i]
      }
    }
  } else if (grepl("cnvkit", as.character(list_of_output_files[1]), fixed = TRUE)) {
    named_files <- vector("list")
    # check for substrings in each file to tell us what the file is
    for (i in 1:length(list_of_output_files)) {
      if (endsWith(list_of_output_files[i], "call.cns") & !endsWith(list_of_output_files[i], "somatic.call.cns")) {
        named_files[["calls"]] <- list_of_output_files[i]
      }
    }
  }
  return(named_files)
}


get_sarek_vcf_file <- function(spn,
                               sampleID,
                               coverage,
                               purity,
                               variant_caller,
                               type,
                               basedir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"){
  file_list <- parse_sarek_variant_called_files(
    get_sarek_variant_called_files(spn = spn,
                                   sampleID = sampleID,
                                   coverage = coverage,
                                   purity = purity,
                                   variant_caller = variant_caller,
                                   type = type,
                                   basedir = basedir)
  )
  if (variant_caller == "strelka") {
    if (type == "tumour") {
      output <- file_list[grep(paste(c("snvs_vcf", "indels_vcf"), collapse = "|"),
                               names(file_list),
                               ignore.case = TRUE)
                          ]
      names(output) <- c("SNV", "INDEL")
    } else if (type == "normal") {
      output <- file_list$variants_vcf
    }
  } else if (variant_caller == "mutect2") {
    output <- file_list$vcf
  } else if (variant_caller == "freebayes") {
    output <- file_list$vcf
  } else if (variant_caller == "haplotypecaller") {
    output <- file_list$vcf
  } else {
    stop("Error: Invalid variant_caller name supplied")
  }
  return(output)
}

get_sarek_cna_file <- function(spn,
                               sampleID,
                               coverage,
                               purity,
                               caller,
                               type=NULL,
                               basedir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/") {
  file_list <- parse_sarek_variant_called_files(
    get_sarek_variant_called_files(spn = spn,
                                   sampleID = sampleID,
                                   coverage = coverage,
                                   purity = purity,
                                   variant_caller = caller,
                                   type = type,
                                   basedir = basedir)
  )
  if (caller == "ascat") {
    output <- file_list
  } else if (caller == "cnvkit") {
    output <- file_list$calls
  } else {
    stop("Error: Invalid caller name supplied")
  }
  return(output)
}

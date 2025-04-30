rm(list = ls())
library(ProCESS)
library(dplyr)
library(bench)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10) {
  stop(paste("Syntax error: ProCESS_seq.R",
	         "<phylo_forest> <SPN> <lot_name>",
	         "<node_local_dir> <output_dir>",
	         "<coverage> <type> <num_of_cores>",
	         "<seed> <purity>"),
       call. = FALSE)
}

phylo_forest_filename <- args[1]
spn_name <- args[2]
lot_name <- args[3]
node_local_dir <- args[4]
output_dir <- args[5]
coverage <- as.double(args[6])
type <- args[7]
num_of_cores <- strtoi(args[8])
seed <- strtoi(args[9])
purity <- as.double(args[10])

if (type == "tumour") {
    seq_tumour <- TRUE
} else if (type == "normal") {
    seq_tumour <- FALSE
    with_preneoplastic <- FALSE
} else if (type == "normal_with_preneoplastic") {
    seq_tumour <- FALSE
    with_preneoplastic <- TRUE
} else {
    stop(paste("The paramter <type> must be among",
               "\"tumour\", \"normal\", and
               \"normal_with_preneoplastic\"."),
         call. = FALSE)
}

merge_sams <- function(output_local_dir, BAM_file,
                       sam_filename_prefix, chromosomes, num_of_cores, resources_dir) {
    
    SAM_files <- ""
    for (i in 1:length(chromosomes)) {
        chr_SAM_file <- file.path(output_local_dir,
                                  paste0(sam_filename_prefix,
                                         chromosomes[i], ".sam"))

        SAM_files <- paste(SAM_files, chr_SAM_file)
    }
    
    name <- paste0(resources_dir, "/out_samtools_merge_", purity, "_", sam_filename_prefix, chromosomes[i])
    cmd <- paste("/usr/bin/time -o", name, "-p -v samtools merge -fc -@", num_of_cores,
		 "-o", BAM_file, SAM_files)
    print(cmd)
    invisible(system(cmd, intern = TRUE))
}

delete_sams <- function(output_local_dir, sam_filename_prefix, chromosomes) {
    for (i in 1:length(chromosomes)) {
        chr_SAM_file <- file.path(output_local_dir,
                                  paste0(sam_filename_prefix,
                                         chromosomes[i], ".sam"))

        unlink(chr_SAM_file)
    }
}

simulate_seq_resources <- function(c,tumour,ref_path,coverage,purity){
  p_info <- ps::ps_handle()
  start_time <- Sys.time()
  initial_cpu <- ps::ps_cpu_times(p_info)
  initial_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
  if (tumour){
    seq_res <- simulate_seq(phylo_forest, reference_genome = ref_path,
                            chromosomes = c,
                            coverage = coverage,
                            purity = purity, 
                            write_SAM = TRUE, 
                            read_size = 150,
                            sequencer = basic_seq,
                            insert_size_mean = 350,
                            insert_size_stddev = 10,
                            output_dir = output_local_dir,
                            include_non_sequenced_mutations = TRUE,
                            update_SAM = TRUE,
                            filename_prefix = sam_filename_prefix,
                            template_name_prefix = paste0(lot_name,'r'),
                            with_normal_sample = FALSE)
  } else{
    cat("simulate_normal_seq")
    seq_res <- simulate_normal_seq(phylo_forest, reference_genome = ref_path,
                                   chromosomes = c,
                                   coverage = coverage,
                                   write_SAM = TRUE, 
                                   read_size = 150,
                                   sequencer = basic_seq,
                                   insert_size_mean = 350,
                                   include_non_sequenced_mutations = TRUE,
                                   insert_size_stddev = 10,
                                   filename_prefix = sam_filename_prefix,
                                   template_name_prefix = paste0(lot_name,'r'),
                                   output_dir = output_local_dir,
                                   with_preneoplastic = with_preneoplastic,
                                   update_SAM = TRUE)
  }

  end_time <- Sys.time()
  final_cpu <- ps::ps_cpu_times(p_info)
  final_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
  elapsed_time <- end_time - start_time
  elapsed_time <- as.numeric(elapsed_time, units = "mins")
  cpu_used <- (final_cpu["user"] + final_cpu["system"]) - (initial_cpu["user"]+ initial_cpu["system"])
  mem_used <- final_mem - initial_mem
  resource_usage <- data.frame(
    elapsed_time_mins =  elapsed_time,
    cpu_time_secs = cpu_used,
    memory_used_MB = mem_used,
    chr = c,
    coverage = coverage,
    purity = purity,
    tumour = tumour
  )
  rownames(resource_usage) <- NULL
  out <- list(mutations=seq_res$mutations,parameters=seq_res$parameters,
              resource_usage = resource_usage)
  return(out)
}

if (!file.exists(node_local_dir)) {
  dir.create(node_local_dir)
}

output_local_dir <- file.path(node_local_dir,
                              paste0(spn_name, "_",
                                     lot_name))
if (file.exists(output_local_dir)) {
  unlink(output_local_dir, recursive=TRUE)
}
dir.create(output_local_dir)

if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

bam_dir <- file.path(output_dir, "BAM")
if (!file.exists(bam_dir)) {
  dir.create(bam_dir)
}

fastq_dir <- file.path(output_dir, "FASTQ")
if (!file.exists(fastq_dir)) {
  dir.create(fastq_dir)
}

resources_dir <- file.path(output_dir, "TIME")
if (!file.exists(resources_dir)) {
  dir.create(resources_dir)
}

data_dir_muts <- file.path(output_dir, "data/mutations")
data_dir_params <- file.path(output_dir, "data/parameters")
data_dir_resources <- file.path(output_dir, "data/resources")

if (!file.exists(data_dir_muts)) {
  dir.create(data_dir_muts,recursive = T)
}

if (!file.exists(data_dir_params)) {
  dir.create(data_dir_params,recursive = T)
}

if (!file.exists(data_dir_resources)) {
  dir.create(data_dir_resources,recursive = T)
}


set.seed(seed)

filename_prefix <- lot_name

sam_filename_prefix <- paste0(filename_prefix, "_chr_")

BAM_filename <- paste0(filename_prefix, ".bam")

BAM_file <- file.path(bam_dir, BAM_filename)
BAM_local_file <- file.path(output_local_dir, BAM_filename)

BAM_done_filename <- file.path(output_dir, paste0(lot_name, "_BAM.done"))

step <- 1

if (!file.exists(BAM_done_filename) || !file.exists(BAM_file)) {
    unlink(BAM_done_filename)

    cat("1. Reading phylogenetic forest...\n")
    phylo_forest <- load_phylogenetic_forest(phylo_forest_filename)
    
    cat("done\n2. Copying reference genome...")
    
    ref_path <- file.path(output_local_dir, "reference.fasta")
    
    invisible(file.copy(phylo_forest$get_reference_path(), ref_path))
    
    cat("done\n3. Simulating reads...\n")
    
    # Simulate sequencing ####
    basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose
    chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr
    cat("parallel::mclapply(",chromosomes,",simulate_seq_resources,mc.cores = 4,tumour = ",
    seq_tumour,",ref_path=",ref_path,",coverage = ",coverage,",purity=",purity,")")
    
    seq_results <- parallel::mclapply(chromosomes, simulate_seq_resources,
                          mc.cores = 4,
                          tumour = seq_tumour, ref_path= ref_path,
                          coverage = coverage, purity = purity)
                          
    seq_results_muts_final <- lapply(1:length(seq_results), function(i) {
        s <- seq_results[[i]]$mutations
        }) %>% do.call(bind_rows, .)


    seq_results_params_final <- lapply(1:length(chromosomes), function(i) {
        pp <- seq_results[[i]]$parameters
        dplyr::tibble(chr = chromosomes[i], parameters = list(pp))
        }) %>% do.call("bind_rows", .)
        
    seq_results_resources_final <- lapply(1:length(seq_results), function(i) {
        s <- seq_results[[i]]$resource_usage
        }) %>% do.call(bind_rows, .)

    saveRDS(seq_results_muts_final,
            file.path(data_dir_muts,
                      paste0("/seq_results_muts_", spn_name,
    			  "_", lot_name, ".rds")))
    saveRDS(seq_results_params_final,
            file.path(data_dir_params,
                      paste0("/seq_results_params_", spn_name,
                  "_", lot_name, ".rds")))
                  
    saveRDS(seq_results_resources_final,
            file.path(data_dir_resources,
                      paste0("/seq_results_resources_", spn_name,
                  "_", lot_name, ".rds")))
    

    cat("done\n4. Building overall BAM file...")
    merge_sams(output_local_dir, BAM_local_file,
               sam_filename_prefix, chromosomes,
    		   num_of_cores, resources_dir)
    
    cat("done\n5. Deleting SAM files...")
    delete_sams(output_local_dir, sam_filename_prefix, chromosomes)
    
    cat("done\n6. Moving the BAM file to output directory...")
    cmd <- paste0("cp ", BAM_local_file, " ", bam_dir, "/")
    invisible(system(cmd, intern = TRUE))

    invisible(file.create(BAM_done_filename))
    cat("done\n")
    remove_local_bam <- TRUE
    step <- 7

} else {
    
    BAM_local_file <- BAM_file
    cat("Found the lot BAM file\n")
    remove_local_bam <- FALSE
    step <- 1
}

cat(paste0(step, ". Splitting BAM file by sample..."))
step <- step + 1

split_bam_by_samples <- function(output_local_dir, BAM_file, remove_local_bam, num_of_cores, resources_dir) {
    name <- strsplit(BAM_file, '/') %>% unlist()
    name <- name[length(name)]
    cmd <- paste0("/usr/bin/time -o ", resources_dir, "/out_samtools_split_", purity,"_",
                  name, " -p -v samtools split -f \"",
                  file.path(output_local_dir,"%*_%!.bam"),
                  "\" ", BAM_file, " -@ ", num_of_cores)
    invisible(system(cmd, intern = TRUE))

    if (remove_local_bam) {
        file.remove(BAM_file)
    }
}
invisible(split_bam_by_samples(output_local_dir, BAM_local_file, remove_local_bam, num_of_cores, resources_dir))
## unlink(bam_dir, recursive = TRUE)

cat(paste0("done\n", step,
           ". Generating the FASTQs and deleting the BAMs..."))
step <- step + 1

BAM_files <- list.files(output_local_dir, pattern = "\\.bam$")

generate_fastq <- function(orig_file, fastq_dir, resources_dir) {
  base_orig_file <- tools::file_path_sans_ext(basename(orig_file))

  file_prefix <- file.path(fastq_dir, base_orig_file)
  R1 <- paste0(file_prefix, ".R1.fastq.gz")
  R2 <- paste0(file_prefix, ".R2.fastq.gz")
  unpaired <- paste0(file_prefix, ".unpaired.fastq.gz")
  singleton <- paste0(file_prefix, ".singleton.fastq.gz")

  prefix <- strsplit(file_prefix,'/') %>% unlist()
  prefix <- prefix[length(prefix)]
  name <- paste0(resources_dir, "/out_fastq_", purity,"_",prefix) 
  cmd <- paste("/usr/bin/time -o", name, "-p -v samtools fastq -@ 20 -c 9 -N -1", R1, "-2", R2, "-0", unpaired, 
               "-s", singleton, orig_file)
  invisible(system(cmd, intern = TRUE))
}

result <- parallel::mclapply(BAM_files, function(c) {
    curr_BAM_file <- file.path(output_local_dir, c)
    if (BAM_file != curr_BAM_file) {
        generate_fastq(curr_BAM_file, output_local_dir, resources_dir)

        unlink(curr_BAM_file)
    }
}, mc.cores = num_of_cores)

cat(paste0("done\n", step,
           ". Moving the FASTQ files to output directory..."))
step <- step + 1

cmd <- paste0("mv ", file.path(output_local_dir, "*.fastq.gz"),
              " ", fastq_dir, "/")
invisible(system(cmd, intern = TRUE))

cat(paste0("done\n", step, ". Removing local files..."))
step <- step + 1

unlink(output_local_dir, recursive = TRUE)

done_filename <- file.path(output_dir, paste0(lot_name, "_final.done"))
invisible(file.create(done_filename))

cat("done\n")

rm(list = ls())
library(rRACES)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10) {
  stop(paste("Syntax error: rRACES_seq.R",
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

merge_sams <- function(output_local_dir, filename_prefix,
                       sam_filename_prefix, chromosomes,
		       num_of_cores) {
    BAM_file <- file.path(output_local_dir,
                          paste0(filename_prefix, ".bam"))
    
    SAM_files <- ""
    for (i in 1:length(chromosomes)) {
        chr_SAM_file <- file.path(output_local_dir,
                                  paste0(sam_filename_prefix,
                                         chromosomes[i], ".sam"))

        SAM_files <- paste(SAM_files, chr_SAM_file)
    }
    
    cmd <- paste("samtools merge -fc -@", num_of_cores,
		 "-o", BAM_file, SAM_files)

    invisible(system(cmd, intern = TRUE))

    return(BAM_file)
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

data_dir <- file.path(output_dir, "data")
if (!file.exists(data_dir)) {
  dir.create(data_dir)
}

set.seed(seed)

filename_prefix <- lot_name

sam_filename_prefix <- paste0(filename_prefix, "_chr_")

cat("1. Reading phylogenetic forest...\n")
phylo_forest <- load_phylogenetic_forest(phylo_forest_filename)

cat("   done\n2. Copying reference genome...")

ref_path <- file.path(output_local_dir, "reference.fasta")

invisible(file.copy(phylo_forest$get_reference_path(), ref_path))

cat("   done\n3. Simulating reads...\n")

# Simulate sequencing ####
#no_error_seq <- ErrorlessIlluminaSequencer()
basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose
chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr
if (seq_tumour) {
  seq_results <- parallel::mclapply(chromosomes, function(c) {
    simulate_seq(phylo_forest, reference_genome = ref_path,
	             chromosomes = c,
                 coverage = coverage,
                 purity = purity, 
                 write_SAM = TRUE, read_size = 150,
                 sequencer = basic_seq,
                 insert_size_mean = 350,
                 insert_size_stddev = 10,
                 output_dir = output_local_dir,
                 update_SAM = TRUE,
                 filename_prefix = sam_filename_prefix,
                 template_name_prefix = paste0(lot_name,'r'),
                 with_normal_sample = TRUE)
  }, mc.cores = num_of_cores)
} else {
  seq_results <- parallel::mclapply(chromosomes, function(c) {
    simulate_normal_seq(phylo_forest, reference_genome = ref_path,
                        chromosomes = c,
                        coverage = coverage,
                        write_SAM = TRUE, read_size = 150,
                        sequencer = basic_seq,
                        insert_size_mean = 350,
                        insert_size_stddev = 10,
                        filename_prefix = sam_filename_prefix,
                 	    template_name_prefix = paste0(lot_name,'r'),
                        output_dir = output_local_dir,
                        with_preneoplastic = with_preneoplastic,
                        update_SAM = TRUE)
  }, mc.cores = num_of_cores)
}
seq_results_final<- do.call("bind_rows", seq_results)
saveRDS(seq_results_final,
        file.path(data_dir,
                  paste0("seq_results_", spn_name,
			  "_", lot_name, ".rds")))

cat("   done\n4. Building overall BAM file...")
BAM_file <- merge_sams(output_local_dir, filename_prefix,
                       sam_filename_prefix, chromosomes,
		               num_of_cores)

split_bam_by_samples <- function(output_local_dir, BAM_file) {
    cmd <- paste0("samtools split -f \"",
                  file.path(output_local_dir,"%*_%!.bam"),
                  "\" ", BAM_file, " -@ ", num_of_cores)
    invisible(system(cmd, intern = TRUE))

    file.remove(BAM_file)
}

cat("done\n5. Splitting BAM file by sample...")
invisible(split_bam_by_samples(output_local_dir, BAM_file))

cat("done\n6. Moving BAM files to output directory...")
cmd <- paste0("mv ", file.path(output_local_dir, "*.bam"),
              " ", bam_dir, "/")
invisible(system(cmd, intern = TRUE))

cat("done\n7. Removing local files...")
unlink(output_local_dir, recursive = TRUE)
cat("done\n")

done_filename <- file.path(output_dir, paste0(lot_name, ".done"))
invisible(file.create(done_filename))

#!/usr/bin/env python3

import os
import sys
import glob
import math
import time
import subprocess
import argparse
import itertools

coverages = [50, 100, 150, 200]
purities = [0.3, 0.6, 0.9]

diluite_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=6:00:00
#SBATCH --mem=20GB

module load R
module load samtools

echo "Rscript diluite_cohort.R ${SPN} ${OUTPUT_DIR} ${FASTQ_DIR} ${TUMOUR_DIR} ${NORMAL_DIR} ${PHYLO_FOREST} ${COVERAGE} ${PURITY} ${ID} 20"

Rscript diluite_cohort.R ${SPN} ${OUTPUT_DIR} ${FASTQ_DIR} ${TUMOUR_DIR} ${NORMAL_DIR} ${PHYLO_FOREST} ${COVERAGE} ${PURITY} ${ID} 20
"""

diluite_R_script="""rm(list = ls())
library(rRACES)

extract_from <- function(orig_file, sampling_fraction, dest_file, seed = 0,
                         sample_name = "") {
  if (sample_name != "") {
    sample_param <- paste("-r", sample_name)
  } else {
    sample_param <- ""
  }
  cmd <- paste("samtools view", sample_param, "--subsample",
               sampling_fraction, "--subsample-seed", seed,
               "-b -o", dest_file, orig_file)

  #print(cmd)
  invisible(system(cmd, intern = TRUE))
}

get_lot_data <- function(BAM_file, sample_name = "") {
  if (sample_name != "") {
    sample_param <- paste("-I", sample_name)
  } else {
    sample_param <- ""
  }
  cmd <- paste("samtools stats", sample_param, BAM_file,
               " | grep '^SN' | cut -f 2-")

  data <- system(cmd, intern = TRUE)

  read_length <- 0
  for (i in seq_len(length(data))) {
    if (startsWith(data[i], "reads paired")) {
      num_of_reads <- as.double(strsplit(data[i], "\t")[[1]][2])
    }

    if (startsWith(data[i], "maximum length")) {
      read_length <- strtoi(strsplit(data[i], "\t")[[1]][2])
    }
  }

  sample_data <- list(num_of_reads = num_of_reads,
                      read_length = read_length)

  return(sample_data)
}

get_lot_number <- function(lot_filename) {
  return(substring(strsplit(lot_filename, split = "_")[[1]][1], 2))
}

read_sample_data <- function(directory, lot_prefix, sample_name) {
  pattern <- paste0(lot_prefix, "0+_", sample_name, "\\\\.bam$")
  lot0 <- list.files(directory, pattern = pattern)

  if (length(lot0) == 0) {
    stop(paste0("File \\"", file.path(directory, pattern),
                "\\" does not exist"), call. = FALSE)
  }

  lot_data <- get_lot_data(file.path(directory, lot0),
                           sample_name = sample_name)

  pattern <- paste0(lot_prefix, "\\\\d+_", sample_name, "\\\\.bam$")

  lot_data <- c(lot_data, lot_files = 1)

  lot_data$lot_files <- list.files(directory, pattern = pattern)

  total_num_of_reads <- (lot_data$num_of_reads *
                         length(lot_data$lot_files))

  lot_data <- c(lot_data,
                total_num_of_reads = total_num_of_reads)

  return(lot_data)
}

generate_fastq <- function(orig_file, fastq_dir) {
  base_orig_file <- tools::file_path_sans_ext(basename(orig_file))

  file_prefix <- file.path(fastq_dir, base_orig_file)
  R1 <- paste0(file_prefix, ".R1.fastq.gz")
  R2 <- paste0(file_prefix, ".R2.fastq.gz")
  unpaired <- paste0(file_prefix, ".unpaired.fastq.gz")
  singleton <- paste0(file_prefix, ".singleton.fastq.gz")

  cmd <- paste("samtools fastq -c 9 -N -1", R1, "-2", R2, "-0", unpaired, 
               "-s", singleton, orig_file)
  invisible(system(cmd, intern = TRUE))
}

write_sarek_line <- function(fastq_dir, orig_filename, SPN, coverage,
                             purity, sample_name, gender, line_counter) {
  base_orig_file <- tools::file_path_sans_ext(basename(orig_filename))
  
  file_prefix <- file.path(fastq_dir, base_orig_file)
  R1 <- paste0(file_prefix, ".R1.fastq.gz")
  R2 <- paste0(file_prefix, ".R2.fastq.gz")

  line <- paste0("L", formatC(line_counter, width = 3,
                              format = "d", flag = "0"))
  cat(paste0("\n", SPN, ",", gender, ",1,", coverage, "X_",
             purity, "p_", sample_name, ",",
             line, ",", R1, ",", R2))
}

create_cohort_dir <- function(output_dir, coverage, purity) {
  if (!file.exists(output_dir)) {
    dir.create(output_dir)
  }

  cohort_dir <- file.path(output_dir, paste0(coverage, "X_",
                                             purity, "p"))

  if (!file.exists(cohort_dir)) {
    dir.create(cohort_dir)
  }

  return(cohort_dir)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10) {
  stop(paste("Syntax error: diluite_cohort.R <SPN>",
             "<output_directory> <FASTQ_dir> <tumour_BAM_directory>",
             "<normal_BAM_directory> <rRACES_phylo_forest>",
             "<coverage> <purity> <id> <num_of_threads>"),
       call. = FALSE)
}

SPN <- args[1]
output_dir <- args[2]
fastq_dir <- args[3]
tumour_dir <- args[4]
normal_dir <- args[5]
forest_file <- args[6]
coverage <- as.double(args[7])
purity <- as.double(args[8])
seed <- args[9]
num_of_threads <- args[10]
build_sarek_file <- TRUE

set.seed(seed)

forest <- load_phylogenetic_forest(forest_file)

genome_len <- max(forest$get_absolute_chromosome_positions()$to)

samples_info <- forest$get_samples_info()

cat(paste0("Handling coverage ", coverage, " and purity ",
           purity, "...\\n Reading normal data..."))
normal_data <- read_sample_data(normal_dir, "s", "normal_sample")
cat("done")

gender <- forest$get_germline_subject()$gender
if (gender == "male") { 
    gender <- "XY"
} else if (gender == "female") { 
    gender <- "XX"
} else {
    stop(paste0("Unsupported germline subject gender ",
                "\\"",gender,"\\"."),
         call. = FALSE)
}

if (!file.exists(fastq_dir)) {
  dir.create(fastq_dir)
}

for (s in seq_len(nrow(samples_info))) {
  sample_info <- samples_info[s, ]

  sample_name <- sample_info$name[[1]]

  cat(paste0("\\n Processing sample \\"", sample_name, "\\""))

  tumour_cells <- sample_info$tumour_cells[[1]]
  equivalent_normal_cells <- sample_info$equivalent_normal_cells[[1]]

  cat(".")
  tumour_data <- read_sample_data(tumour_dir, "t", sample_name)

  cat(".")
  if (tumour_data$read_length == 0) {
    stop(paste("Unknown sample \\"", sample_name, "\\""), call. = FALSE)
  }

  if (purity > 1 || purity < 0) {
    stop("The purity must be a real number in [0,1].", call. = FALSE)
  }

  normal_cells <- tumour_cells * (1 - purity) / purity
  cat(paste("\\n   sample tumour equivalent normal cells:",
            equivalent_normal_cells))
  cat(paste("\\n   normal cells:", normal_cells, "tumour cells:",
            tumour_cells))
  normal_prob <- normal_cells / (equivalent_normal_cells + normal_cells)
  cat(paste("\\n   normal prob:", normal_prob,
            "tumour prob:", (1 - normal_prob)))

  num_of_reads <- round(coverage * genome_len / tumour_data$read_length)

  num_of_normal <- rbinom(1, num_of_reads, normal_prob)
  num_of_tumour <- num_of_reads - num_of_normal

  cat(paste0("\\n   # of reads: ", num_of_reads,
             " (normal:", num_of_normal, " tumour:", num_of_tumour, ")"))

  if (num_of_tumour > tumour_data$total_num_of_reads) {
    stop(paste("Not enough tumour reads for coverage", coverage,
               "and purity", purity, "."), call. = FALSE)
  }

  if (num_of_normal > normal_data$total_num_of_reads) {
    stop(paste("Not enough normal reads for coverage", coverage,
               "and purity", purity, "."), call. = FALSE)
  }

  cohort_dir <- create_cohort_dir(output_dir, coverage, purity)

  cohort_sample_dir <- file.path(cohort_dir, sample_name)
  if (!file.exists(cohort_sample_dir)) {
    dir.create(cohort_sample_dir)
  }

  tumour_extract_done <- file.path(cohort_sample_dir,
                                   paste0("tumour_extract.done"))
  if (!file.exists(tumour_extract_done)) {
    cat("\\n   Extracting tumour reads...")
    tumour_fraction <- num_of_tumour / (tumour_data$total_num_of_reads)

    parallel::mclapply(tumour_data$lot_files, function(orig_file) {
        dest_file <- file.path(cohort_sample_dir, orig_file)
        orig_file <- file.path(tumour_dir, orig_file)
        extract_from(orig_file, tumour_fraction, dest_file,
                     seed = seed, sample_name = sample_name)
    }, mc.cores = num_of_threads)

    invisible(file.create(tumour_extract_done))
    cat("done")
  }

  normal_extract_done <- file.path(cohort_sample_dir,
                                   paste0("normal_extract.done"))
  if (!file.exists(normal_extract_done)) {
    cat("\\n   Extracting normal reads...")
    normal_fraction <- num_of_normal / (normal_data$total_num_of_reads)

    parallel::mclapply(normal_data$lot_files, function(orig_file) {
        dest_file <- file.path(cohort_sample_dir, orig_file)
        orig_file <- file.path(normal_dir, orig_file)

        extract_from(orig_file, normal_fraction, dest_file,
                     seed = seed, sample_name = "normal_sample")
    }, mc.cores = num_of_threads)

    invisible(file.create(normal_extract_done))
    cat("done")
  }

  merge_done_filename <- file.path(cohort_sample_dir, paste0("merge.done"))
  if (FALSE && !file.exists(merge_done_filename)) {
    cat("\\n   Merging normal and tumour reads...")

    final_file <- file.path(cohort_dir, paste0(sample_name, ".cram"))
    cmd <- paste("samtools merge -f -c -O cram --reference",
                forest$get_reference_path(),
                "-@", num_of_threads, "-o", final_file, 
                file.path(cohort_sample_dir, "*.bam"))
    invisible(system(cmd, intern = TRUE))

    invisible(file.create(merge_done_filename))

    cat("done")
  }

  cohort_fastq_dir <- file.path(fastq_dir,
                                paste0(coverage, "X_", purity, "p"))
  if (!file.exists(cohort_fastq_dir)) {
    dir.create(cohort_fastq_dir)
  }

  fastq_sample_dir <- file.path(cohort_fastq_dir, sample_name)
  if (!file.exists(fastq_sample_dir)) {
    dir.create(fastq_sample_dir)
  }

  fastq_done_filename <- file.path(cohort_sample_dir, paste0("fastq.done"))
  if (!file.exists(fastq_done_filename)) {
    cat("\\n   Generating fastq files...")

    parallel::mclapply(tumour_data$lot_files, function(orig_file) {
        generate_fastq(file.path(cohort_sample_dir, orig_file), fastq_sample_dir)
    }, mc.cores = num_of_threads)
    parallel::mclapply(normal_data$lot_files, function(orig_file) {
        generate_fastq(file.path(cohort_sample_dir, orig_file), fastq_sample_dir)
    }, mc.cores = num_of_threads)

    invisible(file.create(fastq_done_filename))
    cat("done")
  }

  sarek_done_filename <- file.path(cohort_sample_dir, paste0("sarek.done"))
  if (build_sarek_file && !file.exists(sarek_done_filename)) {
    cat("\\n   Generating sarek file...")

    sarek_filename <- file.path(cohort_sample_dir,
                                paste0("sarek_",sample_name,".csv"))

    if (!file.exists(sarek_filename)) {
      file.create(sarek_filename)
    }
    sarek_file <- sink(sarek_filename)

    cat("patient,sex,status,sample,lane,fastq_1,fastq_2")
    lot_files <- c(tumour_data$lot_files, normal_data$lot_files)
    for (i in seq_len(length(lot_files))) {
        orig_filename <- lot_files[i]
        invisible(write_sarek_line(fastq_sample_dir, orig_filename,
                                   SPN, coverage, purity, sample_name,
                                   gender, i))
    }
    sink()

    invisible(file.create(sarek_done_filename))
    cat("done")
  }


  cat("\\n done\\ndone")
}

done_filename <- file.path(output_dir, paste0("diluition_", seed, ".done"))
invisible(file.create(done_filename))
"""


def remove_old_done_files(output_dir, lot_prefix):
    done_files = glob.glob(f"{output_dir}/diluition_*.done")
    for done_file in done_files:
        os.unlink(done_file)


def get_completed_jobs(output_dir):
    common_prefix = f"{output_dir}/diluition_"
    common_suffix = '.done'
    done_files = glob.glob(f"{common_prefix}*{common_suffix}")
    prefix_len = len(common_prefix)
    suffix_len = len(common_suffix)

    done_ids = list()
    for done_file in done_files:
        done_ids.append(int(done_file[prefix_len:-suffix_len]))

    return done_ids


if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description=('Diluites tumour and normal '
                                                  + 'BAM and produces FASTQ files'))
    parser.add_argument('SPN', type=str, help='The SPN name (e.g., SPN01)')
    parser.add_argument('phylogenetic_forest', type=str,
                        help = ('The phylogenetic forest absoluted path '
                                + 'for the cluster nodes'))
    parser.add_argument('output_dir', type=str,
                        help = ('The absoluted path of the output'
                                + 'directory'))
    parser.add_argument('fastq_dir', type=str,
                        help = ('The absoluted path of the fastq'
                                + 'directory'))
    parser.add_argument('tumour_dir', type=str,
                        help = ('The absoluted path of the tumour'
                                + 'directory'))
    parser.add_argument('normal_dir', type=str,
                        help = ('The absoluted path of the normal'
                                + 'directory'))
    parser.add_argument('-p', '--partition', type=str, required=True,
                        help="The cluster partition")
    parser.add_argument('-A', '--account', type=str,
                        help="The cluster account")
    parser.add_argument('-j', '--parallel_jobs', type=int, default=20,
                        help="The number of parallel jobs")
    parser.add_argument('-x', '--exclude', type=str, default="",
                        help=("A list of nodes to exclude from the "
                              + "computation"))
    
    args = parser.parse_args()

    if args.account is None:
        process = subprocess.Popen(['whoami'],
                                stdout=subprocess.PIPE)
        account = process.communicate()
    else:
        account = args.account

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    log_dir = '{}/log/'.format(args.output_dir)
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    with open('diluite_cohort.R', 'w') as outstream:
        outstream.write(diluite_R_script)

    with open('diluite_cohort.sh', 'w') as outstream:
        outstream.write(diluite_shell_script)

    params = list(itertools.product(coverages, purities))
    lot_ids = set(range(len(params)))

    completed_ids = get_completed_jobs(args.output_dir)
    submitted = list(completed_ids)
    lot_ids = list(lot_ids.difference(set(completed_ids)))

    zeros = math.ceil(math.log10(len(params)))
    while len(lot_ids) != 0:
        completed_ids = get_completed_jobs(args.output_dir)

        to_be_submitted = (args.parallel_jobs
                           + len(completed_ids)
                           - len(submitted))

        for i in lot_ids[:to_be_submitted]:
            (coverage, purity) = params[i]

            lot_name = 'diluition_{}'.format(i)
            sys.stdout.write('Submitting lot {}...'.format(lot_name))
            sys.stdout.flush()

            cmd = ['sbatch', '--account={}'.format(account),
                   '--partition={}'.format(args.partition),
                   '--job-name={}'.format(lot_name),
                   ('--export=OUTPUT_DIR={},FASTQ_DIR={},TUMOUR_DIR={},'
                    + 'NORMAL_DIR={},PHYLO_FOREST={},COVERAGE={},'
                    + 'PURITY={},SPN={},'
                    + 'ID={}').format(args.output_dir, args.fastq_dir,
                                      args.tumour_dir, args.normal_dir,
                                      args.phylogenetic_forest,
                                      coverage, purity, args.SPN, i),
                   '--output={}/{}.log'.format(log_dir, lot_name),
                   './diluite_cohort.sh']
            if args.exclude != "":
                cmd.insert(-1,"--exclude={}".format(args.exclude))

            subprocess.run(cmd)

            sys.stdout.write('done\n')
            sys.stdout.flush()

        if to_be_submitted>0:
            submitted.extend(lot_ids[:to_be_submitted])
            lot_ids = lot_ids[to_be_submitted:]
        time.sleep(60)  # wait 1 minute

    completed_ids = get_completed_jobs(args.output_dir)
    while (len(completed_ids) != len(submitted)):
        time.sleep(60)
        completed_ids = get_completed_jobs(args.output_dir)

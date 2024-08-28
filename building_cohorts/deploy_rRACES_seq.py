#!/usr/bin/python3

import os
import sys
import math
import glob
import time
import subprocess
import argparse


shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem=80GB

module load R
module load samtools

echo "Rscript rRACES_seq.R ${PHYLO_FOREST} ${SPN} ${LOT} ${NODE_SCRATCH} ${DEST} ${COVERAGE} ${TYPE} 4 ${SEED}"

Rscript rRACES_seq.R ${PHYLO_FOREST} ${SPN} ${LOT} ${NODE_SCRATCH} ${DEST} ${COVERAGE} ${TYPE} 4 ${SEED}

rm -rf ${NODE_SCRATCH}/${SPN}_${LOT}
"""

R_script="""rm(list = ls())
library(rRACES)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  stop(paste("Syntax error: rRACES_seq.R",
	         "<phylo_forest> <SPN> <lot_name>",
	         "<node_local_dir> <output_dir>",
	         "<coverage> <type> <num_of_cores>",
	         "<seed>"),
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
               "\\"tumour\\", \\"normal\\", and
               \\"normal_with_preneoplastic\\"."),
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

cat("1. Reading phylogenetic forest...\\n")
phylo_forest <- load_phylogenetic_forest(phylo_forest_filename)

cat("   done\\n2. Copying reference genome...")

ref_path <- file.path(output_local_dir, "reference.fasta")

invisible(file.copy(phylo_forest$get_reference_path(), ref_path))

cat("   done\\n3. Simulating reads...\\n")

# Simulate sequencing ####
#no_error_seq <- ErrorlessIlluminaSequencer()
basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose
chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr
if (seq_tumour) {
  seq_results <- parallel::mclapply(chromosomes, function(c) {
    simulate_seq(phylo_forest, reference_genome = ref_path,
	             chromosomes = c,
                 coverage = coverage,
                 write_SAM = TRUE, read_size = 150,
                 sequencer = basic_seq,
                 insert_size_mean = 350,
                 insert_size_stddev = 10,
                 output_dir = output_local_dir,
                 update_SAM = TRUE,
                 filename_prefix = sam_filename_prefix,
                 template_name_prefix = paste0(lot_name,'r'),
                 with_normal_sample = FALSE)
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

cat("   done\\n4. Building overall BAM file...")
BAM_file <- merge_sams(output_local_dir, filename_prefix,
                       sam_filename_prefix, chromosomes,
		               num_of_cores)

split_bam_by_samples <- function(output_local_dir, BAM_file) {
    cmd <- paste0("samtools split -f \\"",
                  file.path(output_local_dir,"%*_%!.bam"),
                  "\\" ", BAM_file, " -@ ", num_of_cores)
    invisible(system(cmd, intern = TRUE))

    file.remove(BAM_file)
}

cat("done\\n5. Splitting BAM file by sample...")
invisible(split_bam_by_samples(output_local_dir, BAM_file))

cat("done\\n6. Moving BAM files to output directory...")
cmd <- paste0("mv ", file.path(output_local_dir, "*.bam"),
              " ", bam_dir, "/")
invisible(system(cmd, intern = TRUE))

cat("done\\n7. Removing local files...")
unlink(output_local_dir, recursive = TRUE)
cat("done\\n")

done_filename <- file.path(output_dir, paste0(lot_name, ".done"))
invisible(file.create(done_filename))
"""


def get_lot_prefix(seq_type):
    if seq_type=='normal':
        return 'n'
    if seq_type=='normal_with_preneoplastic':
        return 's'
    if seq_type=='tumour':
        return 't'
    raise TypeError("Only \"tumour\", \"normal\", and "
                    + "\"normal_with_preneoplastic\" are supported")


def get_completed_jobs(output_dir, lot_prefix):
    common_prefix = f"{output_dir}/{lot_prefix}"
    common_suffix = '.done'
    done_files = glob.glob(f"{common_prefix}*{common_suffix}")
    prefix_len = len(common_prefix)
    suffix_len = len(common_suffix)

    done_ids = list()
    for done_file in done_files:
        done_ids.append(int(done_file[prefix_len:-suffix_len]))

    return done_ids


def remove_old_done_files(output_dir, lot_prefix):
    done_files = glob.glob(f"{output_dir}/{lot_prefix}*.done")
    for done_file in done_files:
        os.unlink(done_file)

def get_sample_names(output_dir, first_lot_name):
    common_prefix = f"{output_dir}/{first_lot_name}_"
    common_suffix = ".done"
    BAM_files = glob.glob(f"{common_prefix}*{common_suffix}")
    prefix_len = len(common_prefix)
    suffix_len = len(common_suffix)

    sample_names = list()
    for BAM_file in BAM_files:
        sample_names.append(BAM_file[prefix_len:-suffix_len])

    return sample_names

if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description=('Produces a rRACES sequencing'))
    parser.add_argument('SPN', type=str, help='The SPN name (e.g., SPN01)')
    parser.add_argument('phylogenetic_forest', type=str,
                        help = ('The phylogenetic forest absoluted path '
                                + 'for the cluster nodes'))
    parser.add_argument('output_dir', type=str,
                        help = ('The output directory absoluted path '
                                + 'for the cluster nodes'))
    parser.add_argument('-p', '--partition', type=str, required=True,
                        help="The cluster partition")
    parser.add_argument('-s', '--node_scratch_directory', type=str,
                        default='/local_scratch',
                        help="The nodes' scratch directory")
    parser.add_argument('-A', '--account', type=str,
                        help="The cluster account")
    parser.add_argument('-c', '--coverage', type=float, default=200.0,
                        help="The final sample overall coverage")
    parser.add_argument('-l', '--num_of_lots', type=int, default=50,
                        help="The number of sequencing lots")
    parser.add_argument('-f', '--first_lot_id', type=int, default=0,
                        help="The first lot identifier")
    parser.add_argument('-j', '--parallel_jobs', type=int, default=20,
                        help="The number of parallel jobs")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-n', '--normal', action='store_true',
                       default=False,
                       help=("Sequencing simulation of a normal sample "
                             + "*without* preneoplastic mutations."))
    group.add_argument('-w', '--normal_with_preneoplastic',
                       action='store_true', default=False,
                       help=("Sequencing simulation of a normal sample "
                             + "*with* preneoplastic mutations."))
    parser.add_argument('-x', '--exclude', type=str, default="",
                        help=("A list of nodes to exclude from the "
                              + "computation"))
    parser.add_argument('-F', '--force_completed_jobs', action='store_true',
                        help=("A Boolean flag to force rerun of "
                              + "already completed job."))

    args = parser.parse_args()

    if args.account is None:
        process = subprocess.Popen(['whoami'],
                                stdout=subprocess.PIPE)
        account = process.communicate()
    else:
        account = args.account

    with open('rRACES_seq.R', 'w') as outstream:
        outstream.write(R_script)

    with open('rRACES_seq.sh', 'w') as outstream:
        outstream.write(shell_script)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    log_dir = '{}/log/'.format(args.output_dir)

    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    lot_coverage = args.coverage/args.num_of_lots

    zeros = math.ceil(math.log10(args.num_of_lots))

    if args.normal:
        seq_type = 'normal'
    else:
        if args.normal_with_preneoplastic:
            seq_type = 'normal_with_preneoplastic'
        else:
            seq_type = 'tumour'

    lot_prefix = get_lot_prefix(seq_type)

    if args.force_completed_jobs:
        remove_old_done_files(args.output_dir, lot_prefix)
    
    lot_ids = set(range(args.first_lot_id,
                        args.first_lot_id+args.num_of_lots))

    completed_ids = get_completed_jobs(args.output_dir, lot_prefix)
    submitted = list(completed_ids)
    lot_ids = list(lot_ids.difference(set(completed_ids)))

    while len(lot_ids) != 0:
        completed_ids = get_completed_jobs(args.output_dir, lot_prefix)
        
        to_be_submitted = (args.parallel_jobs
                           + len(completed_ids)
                           - len(submitted))

        for i in lot_ids[:to_be_submitted]:
            lot_name = '{}{}'.format(lot_prefix, str(i).zfill(zeros))

            sys.stdout.write('Submitting lot {}...'.format(lot_name))
            sys.stdout.flush()

            cmd = ['sbatch', '--account={}'.format(account),
                   '--partition={}'.format(args.partition),
                   '--job-name={}_{}'.format(args.SPN, lot_name),
                   ('--export=PHYLO_FOREST={},SPN={},LOT={},DEST={},'
                    + 'COVERAGE={},TYPE={},NODE_SCRATCH={},'
                    + 'SEED={}').format(args.phylogenetic_forest,
                                        args.SPN, lot_name,
                                        args.output_dir, lot_coverage,
                                        seq_type, args.node_scratch_directory,
                                        i),
                   '--output={}/lot_{}.log'.format(log_dir, lot_name),
                   './rRACES_seq.sh']
            if args.exclude != "":
                cmd.insert(-1,"--exclude={}".format(args.exclude))

            subprocess.run(cmd)
            sys.stdout.write('done\n')
            sys.stdout.flush()

        if to_be_submitted>0:
            submitted.extend(lot_ids[:to_be_submitted])
            lot_ids = lot_ids[to_be_submitted:]
        time.sleep(60)  # wait 1 minute

    completed_ids = get_completed_jobs(args.output_dir, lot_prefix)
    while (len(completed_ids) != len(submitted)):
        time.sleep(60)
        completed_ids = get_completed_jobs(args.output_dir, lot_prefix)
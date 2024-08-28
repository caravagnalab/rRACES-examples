#!/usr/bin/python3

import os
import sys
import glob
import time
import argparse
import subprocess
import datetime
import fnmatch

tumour_coverage = 200
mixing_normal_coverage = 200
normal_coverage = 30
num_of_lots = 200
BAM_parallel_jobs = 40
diluition_parallel_jobs = 12

fastq_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=6:00:00
#SBATCH --mem=20GB

module load R
module load samtools

echo "Rscript generate_normal_fastq.R ${BAM_FILE} ${OUTPUT_DIR}"

Rscript generate_normal_fastq.R ${BAM_FILE} ${OUTPUT_DIR}
"""

fastq_R_script="""rm(list = ls())
library(rRACES)

generate_fastq <- function(orig_file, fastq_dir) {
  base_orig_file <- tools::file_path_sans_ext(basename(orig_file))

  file_prefix <- file.path(fastq_dir, base_orig_file)
  R1 <- paste0(file_prefix, ".R1.fastq.gz")
  R2 <- paste0(file_prefix, ".R2.fastq.gz")
  unpaired <- paste0(file_prefix, ".unpaired.fastq.gz")
  singleton <- paste0(file_prefix, ".singleton.fastq.gz")

  cmd <- paste("samtools fastq -@ 20 -c 9 -N -1", R1, "-2", R2, "-0", unpaired, 
               "-s", singleton, orig_file)
  invisible(system(cmd, intern = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(paste("Syntax error: generate_normal_fastq <BAM_file>",
             "<output_dir>"),
       call. = FALSE)
}

bam_file <- args[1]
output_dir <- args[2]

if (!file.exists(output_dir)) {
    dir.create(output_dir)
}

log_dir <- file.path(output_dir, "log")

if (!file.exists(log_dir)) {
    dir.create(log_dir)
}

base_orig_file <- tools::file_path_sans_ext(basename(bam_file))

generate_fastq(bam_file, output_dir)

done_file <- file.path(log_dir, paste0(base_orig_file, ".done"))

invisible(file.create(done_file))
"""


def build_reads(SPN, phylogenetic_forest, output_dir, account, partition,
                num_of_lots, coverage, parallel_jobs, options=[]):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    cmd = (['python3', './deploy_rRACES_seq.py', '-p', f'{partition}',
            '-A', f'{account}', '-l', f'{num_of_lots}', '-c', f'{coverage}',
            '-j', f'{parallel_jobs}', SPN, phylogenetic_forest, output_dir]
            + options)

    subprocess.run(cmd)

    return output_dir


def diluite_reads(SPN, phylogenetic_forest, output_dir, tumour_dir, mixing_normal_dir,
                  account, partition):
    
    diluitions_dir = os.path.join(output_dir, "diluitions")
    if not os.path.exists(diluitions_dir):
        os.mkdir(diluitions_dir)
    
    fastq_dir = os.path.join(output_dir, "fastq")
    if not os.path.exists(fastq_dir):
        os.mkdir(fastq_dir)

    cmd = ['python3', './deploy_diluite_cohort.py', '-p', f'{partition}',
           '-A', f'{account}', '-j', f'{diluition_parallel_jobs}', SPN,
           phylogenetic_forest, diluitions_dir, fastq_dir,
           os.path.join(tumour_dir, "BAM"), os.path.join(mixing_normal_dir, "BAM")]

    subprocess.run(cmd)

    return (diluitions_dir, fastq_dir)


def get_done(directory):
    return glob.glob(os.path.join(directory, "*.done"))


def build_fastq(bam_dir, fastq_dir, account, partition):
    with open('generate_normal_fastq.R', 'w') as outstream:
        outstream.write(fastq_R_script)

    with open('generate_normal_fastq.sh', 'w') as outstream:
        outstream.write(fastq_shell_script)

    if not os.path.exists(fastq_dir):
        os.mkdir(fastq_dir)

    fastq_log_dir = os.path.join(fastq_dir, "log")
    if not os.path.exists(fastq_log_dir):
        os.mkdir(fastq_log_dir)

    bam_files = glob.glob(os.path.join(bam_dir, "*.bam"))

    for bam_file in bam_files:
        basename = os.path.splitext(os.path.basename(bam_file))[0]
        cmd = ['sbatch', f'--account={account}',
                f'--partition={partition}',
                (f'--export=BAM_FILE={bam_file},'
                 + f'OUTPUT_DIR={fastq_dir}'),
                f'--output={fastq_log_dir}/{basename}.log',
                './generate_normal_fastq.sh']
                #'./test.sh']

        subprocess.run(cmd)

    while (len(bam_files) != len(get_done(fastq_log_dir))):
        time.sleep(60)
    
    return fastq_dir

def get_fastq_filenames(bam_dir, fastq_dir):

    normal_fastq_filenames = []
    for bam_file in sorted(glob.glob(os.path.join(bam_dir, "*.bam"))):
        basename = os.path.splitext(os.path.basename(bam_file))[0]
        common_name = os.path.join(fastq_dir, basename)

        normal_fastq_filenames.append((f'{common_name}.R1.fastq.gz',
                                       f'{common_name}.R2.fastq.gz'))

    return normal_fastq_filenames

def get_cohort_in(directory):
    cohorts = set()
    for cohort in fnmatch.filter(os.listdir(directory), "*X_*p"):
        if os.path.isdir(os.path.join(directory, cohort)):
            cohorts.add(cohort)
    
    return cohorts

def build_sarek_files(output_dir, fastq_dir, diluitions_dir, normal_dir):
    sarek_dir = os.path.join(output_dir, "sarek")

    if (not os.path.exists(sarek_dir)):
        os.mkdir(sarek_dir)

    path = os.path.normpath(diluitions_dir)
    diluitions_dir_depth = len(path.split(os.sep))

    normal_fastq_filenames = get_fastq_filenames(os.path.join(normal_dir, "BAM"),
                                                 os.path.join(fastq_dir, "normal"))

    #cohorts = set()
    #for root, dirnames, filenames in os.walk(diluitions_dir):
    #    for filename in fnmatch.filter(filenames, "sarek*.csv"):
    #        path = os.path.normpath(root)
    #        cohorts.add(path.split(os.sep)[diluitions_dir_depth])

    for cohort in get_cohort_in(diluitions_dir):
        with open(os.path.join(sarek_dir, cohort + ".csv"), 'w') as outstream:
            outstream.write('patient,sex,status,sample,'
                            + 'lane,fastq_1,fastq_2\n')

            for root, dirnames, filenames in os.walk(diluitions_dir):
                for filename in fnmatch.filter(filenames, "sarek_*.csv"):
                    with open(os.path.join(root, filename)) as instream:
                        line_num = 0
                        for line in instream:
                            if line_num>0:
                                outstream.write(f'{line.strip()}\n')
                                if (line_num == 1):
                                    (SPN, gender) = line.strip().split(',')[:2]

                            line_num += 1

            line = 1
            for fastq_filenames in normal_fastq_filenames:
                outstream.write(f'{SPN},{gender},0,normal_sample,'
                                + f'L{str(line).zfill(3)},{fastq_filenames[0]},'
                                + f'{fastq_filenames[1]}\n')
                line += 1

def build_tumour_reads(SPN, phylogenetic_forest, output_dir, account, partition):
    return build_reads(SPN, phylogenetic_forest, os.path.join(output_dir, "tumour"),
                       account, partition, tumour_coverage, tumour_coverage, 
                       min(BAM_parallel_jobs, tumour_coverage))

def build_mixing_normal_reads(SPN, phylogenetic_forest, output_dir, account, partition):
    return build_reads(SPN, phylogenetic_forest, os.path.join(output_dir, "mixing_normal"),
                       account, partition, mixing_normal_coverage, mixing_normal_coverage,
                       min(BAM_parallel_jobs, mixing_normal_coverage), options=['-w'])

def build_normal_reads(SPN, phylogenetic_forest, output_dir, account, partition):
    return build_reads(SPN, phylogenetic_forest, os.path.join(output_dir, "normal"),
                       account, partition, normal_coverage, normal_coverage,
                       min(BAM_parallel_jobs, normal_coverage), options=['-n'])


if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description=('Generates the cohorts of a tumour'))
    parser.add_argument('SPN', type=str, help='The SPN name (e.g., SPN01)')
    parser.add_argument('phylogenetic_forest', type=str,
                        help = ('The phylogenetic forest absoluted path '
                                + 'for the cluster nodes'))
    parser.add_argument('output_dir', type=str,
                        help = ('The output directory absoluted path '
                                + 'for the cluster nodes'))
    parser.add_argument('-p', '--partition', type=str, required=True,
                        help="The cluster partition")
    parser.add_argument('-A', '--account', type=str,
                        help="The cluster account")

    args = parser.parse_args()

    if args.account is None:
        process = subprocess.Popen(['whoami'],
                                stdout=subprocess.PIPE)
        account = process.communicate()
    else:
        account = args.account

    output_dir = os.path.abspath(args.output_dir)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    execution_track = []

    execution_track.append(f"Starting at {datetime.datetime.now()}\n")
    phylogenetic_forest = os.path.abspath(args.phylogenetic_forest)

    tumour_dir = build_tumour_reads(args.SPN, phylogenetic_forest,
                                    output_dir, account, args.partition)
    execution_track.append(f"Tumour BAMs built at {datetime.datetime.now()}\n")

    mixing_normal_dir = build_mixing_normal_reads(args.SPN, phylogenetic_forest,
                                                  output_dir, account, args.partition)
    execution_track.append(f"Mixing normal BAMs built at {datetime.datetime.now()}\n")
    
    normal_dir = build_normal_reads(args.SPN, phylogenetic_forest,
                                    output_dir, account, args.partition)
    execution_track.append(f"Normal BAMs built at {datetime.datetime.now()}\n")

    (diluitions_dir, fastq_dir) = diluite_reads(args.SPN, phylogenetic_forest, output_dir,
                                                tumour_dir, mixing_normal_dir, account,
                                                args.partition)
    execution_track.append(f"Cohorts built at {datetime.datetime.now()}\n")

    build_fastq(os.path.join(normal_dir, "BAM"),
                os.path.join(fastq_dir, "normal"), account, args.partition)
    execution_track.append("Normal FASTQs built at "
                           + f"{datetime.datetime.now()}\n")

    build_sarek_files(output_dir, fastq_dir, diluitions_dir, normal_dir)
    execution_track.append("Sarek files built at "
                           + f"{datetime.datetime.now()}\n")

    with open(f'{os.path.join(output_dir, "cohorts.done")}', 'w') as outstream:
        for line in execution_track:
            outstream.write(line)

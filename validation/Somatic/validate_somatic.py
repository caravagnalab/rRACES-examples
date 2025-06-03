#!/usr/bin/python3

import os
import sys
import math
import glob
import time
import subprocess
import argparse

somatic_processing_shell_script="""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
# Load R module
module load R/4.4.1 


echo "Processing chromosome ${CHROMOSOME}"

# Run R script with the current chromosome

Rscript preprocess.R --spn_id ${SPN} --purity ${PURITY} --coverage ${COVERAGE} --chr ${CHROMOSOME}
"""


somatic_prepare_report_shell_script="""#!/bin/bash
#SBATCH --mem=8GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

# Load R module
module load R/4.4.1
Rscript prepare_report.R --spn_id ${SPN} --purity ${PURITY} --coverage ${COVERAGE}
"""

cna_report_shell_script="""#!/bin/bash
#SBATCH --job-name=validate_CNA
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=1024M

module load R/4.4.1
Rscript Validate_CNA_calls.R --spn_id ${SPN} --sample_id ${SAMPLE} --purity ${PURITY} --coverage ${COVERAGE}  
"""


if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description=('Produces the cohorts of a SPN'))
    parser.add_argument('SPN', type=str, help='The SPN name (e.g., SPN01)')
    parser.add_argument('coverage', type=int, help='The coverage')
    parser.add_argument('purity', type=str, help='The purity')
    parser.add_argument('-P', '--partition', type=str, required=True,
                        help="The cluster partition")
    parser.add_argument('-A', '--account', type=str, required=True,
                        help="The cluster account")

    

    args = parser.parse_args()

    if args.account is None:
        process = subprocess.Popen(['whoami'],
                                stdout=subprocess.PIPE)
        account = process.communicate()
    else:
        account = args.account
    
    scout_dir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
    
    with open(os.path.join(scout_dir, args.SPN, "process", "subject_gender.txt")) as gender_file:
    gender = gender_file.read().strip()

    # Determine chromosomes based on gender
    if gender == "XX":
        chromosomes = list(range(1, 23)) + ['X']
    elif gender == "XY":
        chromosomes = list(range(1, 23)) + ['X', 'Y']
    else:
        raise ValueError(f"Unexpected gender value: {gender}")
    
    curr_dir = os.getcwd()
    log_dir = '{}/out/'.format(curr_dir)

    with open('run_processing.sh', 'w') as outstream:
        outstream.write(somatic_processing_shell_script)
    
    for chr in chromosomes:

        cmd = ['sbatch', '--account={}'.format(account),
            '--partition={}'.format(args.partition),
            '--job-name=parsing_{}_{}_{}_{}'.format(args.SPN, args.coverage,args.purity,chr),
            ('--export=SPN={},COVERAGE={},PURITY={},CHROMOSOME={}').format(args.SPN,
                                                args.coverage, args.purity,chr),
            '--output={}/parsing_chr{}.log'.format(log_dir, chr),
            './run_processing.sh'] 

        subprocess.run(cmd)
    
    # Run final somatic report

    with open('run_reports.sh', 'w') as outstream:
        outstream.write(somatic_prepare_report_shell_script)
    cmd = ['sbatch', '--account={}'.format(account),
        '--partition={}'.format(args.partition),
        '--job-name=somatic_report_{}_{}_{}'.format(args.SPN, args.coverage,args.purity),
        ('--export=SPN={},COVERAGE={},PURITY={}').format(args.SPN,
                                            args.coverage, args.purity),
        '--output={}/somatic_report_{}_{}_{}.log'.format(log_dir, args.SPN,args.coverage,args.purity),
        './run_reports.sh']
    subprocess.run(cmd)

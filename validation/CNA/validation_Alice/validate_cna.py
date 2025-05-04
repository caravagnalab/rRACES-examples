#!/usr/bin/python3

import os
import sys
import math
import glob
import time
import subprocess
import argparse
import fnmatch

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

Rscript preprocess.R --spn_id ${SPN} --purity ${PURITY} --coverage ${COVERAGE} --chr chr${CHROMOSOME}
"""


somatic_prepare_report_shell_script="""#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out/report_prep.log
#SBATCH --job-name=report

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

    chromosomes = list(range(1, 23)) + ['X', 'Y']
    curr_dir = os.getcwd()
    log_dir = '{}/out/'.format(curr_dir)

    # with open('run_processing.sh', 'w') as outstream:
    #     outstream.write(somatic_processing_shell_script)
    
    # for chr in chromosomes:

    #     cmd = ['sbatch', '--account={}'.format(account),
    #         '--partition={}'.format(args.partition),
    #         '--job-name=parsing_{}_{}_{}_{}'.format(args.SPN, args.coverage,args.purity,chr),
    #         ('--export=SPN={},COVERAGE={},PURITY={},CHROMOSOME={}').format(args.SPN,
    #                                             args.coverage, args.purity,chr),
    #         '--output={}/parsing_chr{}.log'.format(log_dir, chr),
    #         './run_processing.sh'] 

    #     subprocess.run(cmd)
    
    ## Run final somatic report

    # with open('run_reports.sh', 'w') as outstream:
    #     outstream.write(somatic_prepare_report_shell_script)
    # cmd = ['sbatch', '--account={}'.format(account),
    #     '--partition={}'.format(args.partition),
    #     '--job-name=somatic_report_{}_{}_{}'.format(args.SPN, args.coverage,args.purity),
    #     ('--export=SPN={},COVERAGE={},PURITY={}').format(args.SPN,
    #                                         args.coverage, args.purity),
    #     '--output={}/somatic_report_{}_{}_{}.log'.format(log_dir, args.SPN,args.coverage,args.purity),
    #     './run_reports.sh']
    # subprocess.run(cmd)
    
     ## Run final CNA report

    base_scout_dir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
    cna_dir = os.path.join(base_scout_dir, args.SPN, "process", "cna_data")
    cna_files = [f for f in os.listdir(cna_dir) if fnmatch.fnmatch(f, '*_cna.rds')]
    sample_ids = [f.replace('_cna.rds', '') for f in cna_files]
    print(sample_ids)

    with open('validate_cna.sh', 'w') as outstream:
        outstream.write(cna_report_shell_script)
    
    for sample in sample_ids:
        cmd = ['sbatch', '--account={}'.format(account),
            '--partition={}'.format(args.partition),
            '--job-name=cna_report_{}_{}_{}'.format(sample, args.coverage,args.purity),
            ('--export=SPN={},SAMPLE={},COVERAGE={},PURITY={}').format(args.SPN,sample,
                                                args.coverage, args.purity),
            '--output={}/cna_report_{}_{}_{}.log'.format(log_dir, sample,args.coverage,args.purity),
            './validate_cna.sh']
        subprocess.run(cmd)

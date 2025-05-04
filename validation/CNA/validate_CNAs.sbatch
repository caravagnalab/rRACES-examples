#!/bin/bash
#SBATCH --job-name=validate_CNA
#SBATCH --mail-type=ALL
#SBATCH --partition=EPYC
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=validate_CNA_%j.log
#SBATCH --mem-per-cpu=1024M

module load R
Rscript Validate_CNA_calls.R --sample_id SPN01_1.1 --purity 0.6p --coverage 100x  


Rscript Validate_CNA_calls.R --purity 0.3p --coverage 50x

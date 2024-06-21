#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=sequencing
#SBATCH --nodes=2
#SBATCH --mem=200gb
#SBATCH --time=03:00:00
#SBATCH --output=sequencing_%j.out

cd /u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06

module load R
Rscript sequencing.R

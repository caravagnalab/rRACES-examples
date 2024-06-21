#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=mutations
#SBATCH --nodes=2
#SBATCH --mem=200gb
#SBATCH --time=03:00:00
#SBATCH --output=mutations_%j.out

cd /u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06

module load R
Rscript mutations.R

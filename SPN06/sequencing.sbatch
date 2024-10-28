#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=SPN6_SEQ
#SBATCH --nodes=2
#SBATCH --mem=200gb
#SBATCH --time=100:00:00
#SBATCH --output=/u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06/logs/sequencing_%j.out

cd /u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06

module load R
Rscript sequencing.R

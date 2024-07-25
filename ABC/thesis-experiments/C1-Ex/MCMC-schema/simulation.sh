#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --account=cdslab
#SBATCH --job-name=rR_c1MC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=rRACES_sim_%j.out
#SBATCH --error=error_rRACES_%J

module load R/4.2.3
Rscript script.R

#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out/report_prep.log
#SBATCH --job-name=report

# Load R module
module load R
Rscript prepare_report.R

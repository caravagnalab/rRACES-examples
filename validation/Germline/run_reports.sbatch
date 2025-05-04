#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=report.log
#SBATCH --job-name=report

# Load R module
module load R/4.4.1

SPN='SPN03'

Rscript compare.R -s $SPN -t 'freebayes' 
Rscript compare.R -s $SPN -t 'haplotypecaller'
Rscript compare.R -s $SPN -t 'strelka'
#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out/report_prep.log
#SBATCH --job-name=report

# Load R module
module load R/4.4.1

SPN='SPN03'
COV='100'
PUR='0.6'

Rscript prepare_report.R --spn_id $SPN --purity $PUR --coverage $COV

#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out/report_prep.log
#SBATCH --job-name=report

# Load R module
conda activate process

SPN='SPN04'
COV='50'
PUR='0.9'

Rscript prepare_report.R --spn_id $SPN --purity $PUR --coverage $COV

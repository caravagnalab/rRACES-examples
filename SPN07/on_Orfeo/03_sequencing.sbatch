#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=03_sequencing
#SBATCH --nodes=2
#SBATCH --mem=500gb
#SBATCH --time=72:00:00
#SBATCH --output=03_sequencing_%j.out

cd /u/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo
module load R
Rscript 03_sequencing.R

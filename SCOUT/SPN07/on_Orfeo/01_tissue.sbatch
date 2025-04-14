#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=01_tissue
#SBATCH --nodes=2
#SBATCH --mem=200gb
#SBATCH --time=03:00:00
#SBATCH --output=01_tissue_%j.out

cd /u/cdslab/antonelloa/ProCESS-examples/SPN07/on_Orfeo
module load R
Rscript 01_tissue.R

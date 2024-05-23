#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --account=cdslab
#SBATCH --job-name=01_tissue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=03:00:00
#SBATCH --output=01_tissue_%j.out

cd /u/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo
module load R
Rscript 01_tissue.R

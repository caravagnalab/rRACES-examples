#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --account=cdslab
#SBATCH --job-name=tissue
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=03:00:00
#SBATCH --output=tissue_%j.out

cd /u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06

module load R
Rscript tissue.R

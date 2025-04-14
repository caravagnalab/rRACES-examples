#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=tissue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=10
#SBATCH --mem=30gb
#SBATCH --time=1:00:00
#SBATCH --output=/u/cdslab/ahaghighi/scratch/packages/ProCESS-examples/SPN06/logs/tissue_%j.log

cd /u/cdslab/ahaghighi/scratch/packages/ProCESS-examples/SPN06/

#conda deactivate

module purge
module load R


#export PATH=/u/cdslab/ahaghighi/scratch/anaconda3/bin:$PATH
#eval "$(conda shell.bash hook)"
#conda activate races-env

Rscript tissue.R

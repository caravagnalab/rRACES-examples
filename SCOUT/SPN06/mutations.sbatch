#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=mutations
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=10
#SBATCH --mem=200gb
#SBATCH --time=03:00:00
#SBATCH --output=/u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06/logs/mutations_%j.log

cd /u/cdslab/ahaghighi/scratch/packages/rRACES-examples/SPN06

module purge
module load R
Rscript mutations.R

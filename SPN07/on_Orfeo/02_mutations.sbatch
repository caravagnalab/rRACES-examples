#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --account=cdslab
#SBATCH --job-name=02_mutations
#SBATCH --nodes=2
#SBATCH --mem=200gb
#SBATCH --time=10:00:00
#SBATCH --output=02_mutations_%j.out

#cd /u/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo
#cd /orfeo/cephfs/scratch/cdslab/shared/races/GRCh38
module load R
Rscript /u/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo/02_mutations.R

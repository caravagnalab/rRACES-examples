#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --account=cdslab
#SBATCH --job-name=02_mutations
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=03:00:00
#SBATCH --output=02_mutations_%j.out

cd /u/cdslab/antonelloa/rRACES-examples/SPN07/on_Orfeo
module load R
Rscript 02_mutations.R

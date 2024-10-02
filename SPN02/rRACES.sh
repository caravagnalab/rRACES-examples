#!/bin/bash
#SBATCH --job-name=sequencing_100purity
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=400gb
#SBATCH --time=70:00:00
#SBATCH --account=cdslab
#SBATCH --output=sequencing_highpurity.out
#SBATCH --error=sequencing_highpurity.err

module load R/4.3.3

Rscript rRACES-examples/SPN02/03_sequencing.R

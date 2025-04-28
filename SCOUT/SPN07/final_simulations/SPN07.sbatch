#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 300gb
#SBATCH --time=10:00:00
#SBATCH --output=ProCESS.out
#SBATCH --error=ProCESS.err

module load singularity
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.0.0.sif"
# change with your path to the simulate_tissue.R and simulate_mutation.R scripts
base="/orfeo/scratch/cdslab/antonelloa/ProCESS-examples/SCOUT/SPN07/final_simulations"

# singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_tissue.R
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_mutations.R

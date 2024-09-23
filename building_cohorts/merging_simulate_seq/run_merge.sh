#!/bin/bash
#SBATCH --partition=THIN #check if is it fine to run it on EPYC node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8 #24 when running simulate_seq
#SBATCH --mem=80G #400G when running simulate_seq
#SBATCH --time=72:00:00
#SBATCH -A cdslab
#SBATCH --output=merge_rds_%j.out
#SBATCH --error=merge_rds_%j.err

module load R/4.3.3
l=$1 ## lot range, 1:10-11:21-...
c=$2 ## lot type single if running the first merging, if final run the final merge
dir=$3 ## directory where orignal and chunk rds files are present
chrom=$4 ## chromosome (e.g 22), this will parallelize in terms of chromsomes
printf "Rscript merge_rds.R --lot_range ${l} --lot_type ${c} --rds_dir ${dir} --chrom ${chrom}""\n"
Rscript merge_rds.R --lot_range ${l} --lot_type ${c} --rds_dir ${dir} --chrom ${chrom}

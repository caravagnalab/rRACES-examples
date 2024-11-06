#!/bin/bash
#SBATCH --partition=THIN
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=72:00:00
#SBATCH --output=sam_to_fastq_normal_%A_%a.out
#SBATCH --error=sam_to_fastq_normal_%A_%a.err
#SBATCH --array=1-24

module load samtools

## Set up info
# Set the directory where sam files for the normal are stored
#races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/FINAL_DATA/sequencing_30X_basic_error_paired_100_1normal/"
races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/CHECK_PURITY/sequencing_30X_basic_error_paired_350_1normal_new_1_germline/"
sam_files=($(ls $races_dir/*.sam | grep chr | rev | cut -f 1 -d "/"| rev |cut -f 1 -d "."))
sam_file=${sam_files[$SLURM_ARRAY_TASK_ID-1]} ## split by chromosomes

## Define the array of samples ids present in sam header
samples=(normal_sample)


## SAM_TO_FASTQ ###
#### inside the for loop try to sbatch a sh file with command line params
#### samtools_fastq.sh script will perform the conversion for each chr and each sample
for i in ${samples[@]}; do
	sbatch samtools_fastq.sh ${sam_file} $i $races_dir
done

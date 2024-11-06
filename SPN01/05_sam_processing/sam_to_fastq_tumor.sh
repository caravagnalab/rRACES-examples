#!/bin/bash
#SBATCH --partition=THIN
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16   #1
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --output=sam_to_fastq_tumor_%A_%a.out
#SBATCH --error=sam_to_fastq_tumor_%A_%a.err
#SBATCH --array=1-24

module load samtools

## Set up info
#races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/FINAL_DATA/sequencing_100X_basic_error_paired_350_3tumor_new1"
races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/CHECK_PURITY/sequencing_100X_basic_error_paired_350_1tumor_new_1"
sam_files=($(ls $races_dir/*.sam | grep chr | rev | cut -f 1 -d "/"| rev |cut -f 1 -d "."))
sam_file=${sam_files[$SLURM_ARRAY_TASK_ID-1]} ## split by chromosomes

## Define the array of samples ids present in sam header
#samples=(SPN01_Sample_1 SPN01_Sample_2 SPN01_Sample_3)
#cd $races_dir
#time samtools split -f %\!/%*_%\!.sam --threads 8 ${sam_file}.sam
samples=(SPN01_Sample_1)
cd /u/cdslab/ggandolfi/scratch/prj_races/rRACES-examples/SPN01/05_sam_processing
## SAM_TO_FASTQ ###
for i in ${samples[@]}; do
   sbatch samtools_fastq.sh ${sam_file} $i $races_dir
done
#
#
#
### inside the for loop try to sbatch a sh file with command line params

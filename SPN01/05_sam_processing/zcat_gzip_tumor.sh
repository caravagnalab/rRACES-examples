#!/bin/bash
#SBATCH --partition=THIN #EPYC
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16 #16
#SBATCH --mem=200G #200G
#SBATCH --time=72:00:00
#SBATCH --output=zcat_normal_%A_%a.out
#SBATCH --error=zcat_normal_%A_%a.err
#SBATCH --array=1-3

#module load samtools
#module load gatk
#module load picard

## Set up info
races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/CHECK_PURITY/sequencing_100X_basic_error_paired_350_1tumor_new_1"
#sam_files=($(ls $races_dir | grep chr | cut -f 1 -d "."))
#sam_file=${sam_files[$SLURM_ARRAY_TASK_ID-1]} ## split by chromosomes

## Define the array of samples ids present in sam header

samples=(SPN01_Sample_1)
#samples=(SPN01_Sample_1 SPN01_Sample_2 SPN01_Sample_3)
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}
#rm $races_dir/*paired*
cd $races_dir

time cat ${sample}/chr_*R1.fastq > ${sample}/${sample}_R1.fastq
#rm -f ${sample}/chr_*R1.fastq
time pigz -p 16 ${sample}/${sample}_R1.fastq ## cpus-per-task
time cat ${sample}/chr_*R2.fastq > ${sample}/${sample}_R2.fastq
#rm -f ${sample}/chr_*R2.fastq
time pigz -p 16 ${sample}/${sample}_R2.fastq ## cpus-per-task

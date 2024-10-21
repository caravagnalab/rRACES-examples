#!/bin/bash
#SBATCH --partition=THIN #EPYC
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8 #16
#SBATCH --mem=50G #200G
#SBATCH --time=72:00:00
#SBATCH --output=zcat_normal_%A_%a.out
#SBATCH --error=zcat_normal_%A_%a.err
#SBATCH --array=1

#module load samtools
#module load gatk
#module load picard

## Set up info
#races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/FINAL_DATA/sequencing_30X_basic_error_paired_100_1normal"
races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/CHECK_PURITY/sequencing_30X_basic_error_paired_350_1normal_new_1_germline"
## Define the array of samples ids present in sam header

samples=(normal_sample)
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}
cd $races_dir

time cat chr_*${sample}/*R1.fastq > ${sample}/${sample}_R1.fastq
rm -f chr_*R1.fastq
time pigz -p 16 ${sample}/${sample}_R1.fastq ## cpus-per-task
time cat chr_*${sample}/*R2.fastq > ${sample}/${sample}_R2.fastq
rm -f chr_*R2.fastq
time pigz -p 16 ${sample}/${sample}_R2.fastq ## cpus-per-task

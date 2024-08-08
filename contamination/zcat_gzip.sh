#!/bin/bash
#SBATCH --partition=THIN #EPYC
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16 #16
#SBATCH --mem=200G #200G
#SBATCH --time=72:00:00
#SBATCH --output=zcat_gzip_%A_%a.out
#SBATCH --error=zcat_gzip_%A_%a.err

sample_list="$1"
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $sample_list)

races_dir="$2"
purity="$3"
frac=($(echo $purity | cut -f 2 -d "."))
coverage="$4"
## Define the array of samples ids present in sam header

fastq_dir=$races_dir/downsampling/tmp/${sample}/bam2fastq
time cat $fastq_dir/chr_*R1.fastq > $races_dir/downsampling/purity_$frac\_coverage_$coverage/${sample}/${sample}_R1.fastq
#rm -f ${sample}/chr_*R1.fastq
time pigz -p 16 $races_dir/downsampling/purity_$frac\_coverage_$coverage/${sample}/${sample}_R1.fastq  ## cpus-per-task

time cat $fastq_dir/chr_*R2.fastq > $races_dir/downsampling/purity_$frac\_coverage_$coverage/${sample}/${sample}_R2.fastq
#rm -f ${sample}/chr_*R1.fastq
time pigz -p 16 $races_dir/downsampling/purity_$frac\_coverage_$coverage/${sample}/${sample}_R2.fastq  ## cpus-per-task

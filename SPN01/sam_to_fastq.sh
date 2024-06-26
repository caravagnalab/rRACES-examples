#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --output=sam_to_bam_%A_%a.out
#SBATCH --error=sam_to_bam_%A_%a.err
#SBATCH --array=1-24

module load samtools
#module load gatk
module load picard

## Set up info
#races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/sequencing_80X_with_error_paired"
races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/sequencing_80X_with_error_paired_new"
sam_files=($(ls $races_dir | grep chr | cut -f 1 -d "."))
sam_file=${sam_files[$SLURM_ARRAY_TASK_ID-1]} ## split by chromosomes

## Define the array of samples ids present in sam header

# samples=(SPN01_Sample_1 normal_sample)

cd $races_dir

samtools split -f %\!/%*_%\!.sam --threads 4 ${sam_file}.sam
#

## SAM_TO_FASTQ ###
for i in ${samples[@]}; do
       java -jar $picard SamToFastq I=$i/${sam_file}_$i.sam \
	       FASTQ=$i/${sam_file}_$i.R1.fastq \
	       SECOND_END_FASTQ=$i/${sam_file}_$i.R2.fastq \
	       UNPAIRED_FASTQ=$i/${sam_file}_$i.unpaired.fastq \
	       VALIDATION_STRINGENCY=SILENT
       gzip $i/${sam_file}_$i.R1.fastq
       gzip $i/${sam_file}_$i.R2.fastq
done

## merge all chromsomes into a unique bam file
#samtools merge -c -o $races_dir/splitted_sam/${sample}\_merged_sorted.bam $races_dir/splitted_sam/chr_*${sample}_reheadered_sorted.bam

## rehader 
#samtools reheader -c 'sed 's/SN:/SN:chr/g'' $races_dir/splitted_sam/${sample}\_merged_sorted.bam > $races_dir/splitted_sam/${sample}_merged_sorted_chrom.bam
#samtools reheader -c 'sed 's/normal_sample/SPN01_normal_sample/g'' $races_dir/splitted_sam/${sample}_merged_sorted_chrom.bam >$races_dir/splitted_sam/${sample}_merged_sorted_chrom_1.bam
#samtools index $races_dir/splitted_sam/${sample}_merged_sorted_chrom.bam

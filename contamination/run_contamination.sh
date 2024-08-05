#!/bin/bash
#SBATCH --partition=THIN
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --output=contamination_%A_%a.out
#SBATCH --error=contamination_%A_%a.err
#SBATCH --array=1-24

module load samtools

## Set up info

races_dir="path/to/dir/where/sam/from/races/tumor/are/sequenced/"
races_dir_n="path/to/dir/where/sam/from/races/normal/are/sequenced/"
sam_files=($(ls $races_dir/*.sam | grep chr | rev | cut -f 1 -d "/"| rev |cut -f 1 -d "."))
sam_file=${sam_files[$SLURM_ARRAY_TASK_ID-1]} ## split by chromosomes
purity=0.20
frac=($(echo $purity | cut -f 2 -d "."))
original_depth_tumor=200
original_depth_normal=160

new_depth_tumor=100


samtools sort $races_dir/${sam_file}.sam > $races_dir/${sam_file}.sorted.bam
samtools sort $races_dir_n/${sam_file}.sam >  $races_dir_n/${sam_file}.sorted.bam
samtools addreplacerg -r ID:SPN01_Sample_1_downsampled -r SM:SPN01_Sample_1_downsampled -r PL:ILLUMINA -o $races_dir/${sam_file}.sorted.downsampled.ID.bam $races_dir/${sam_file}.sorted.bam ## required to have the same ID once merged
samtools addreplacerg -r ID:SPN01_Sample_1_downsampled -r SM:SPN01_Sample_1_downsampled -r PL:ILLUMINA -o $races_dir_n/${sam_file}.sorted.downsampled.ID.bam $races_dir_n/${sam_file}.sorted.bam ## required to have the same ID once merged

/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/contamination/downsampling.sh -t $races_dir/${sam_file}.sorted.downsampled.ID.bam -n $races_dir_n/${sam_file}.sorted.downsampled.ID.bam -i $original_depth_tumor -j $original_depth_normal -p $purity -d $new_depth_tumor -o $races_dir/${sam_file}.contaminated.$frac.bam


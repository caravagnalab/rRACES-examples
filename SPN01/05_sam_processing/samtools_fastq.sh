#!/bin/bash
#SBATCH --partition=THIN
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --output=samtools_fastq_%J.out
#SBATCH --error=samtools_fastq_%J.err

module load samtools

races_dir=$3
cd $races_dir
sample_id=$2
chrom=$1
#cd $sample_id
#date
#printf "samtools fastq -N $chrom\_$sample_id.sam -1 $chrom\_$sample_id.R1.fastq -2 $chrom\_$sample_id.R2.fastq""\n"
#time samtools fastq -N $chrom.sam -1 $chrom\_$sample_id.R1.fastq -2 $chrom\_$sample_id.R2.fastq
#time samtools sort -n $chrom\_$sample_id.sam | \
#	samtools fastq -1 $chrom\_$sample_id.R1.fastq -2 $chrom\_$sample_id.R2.fastq -0 $chrom\_$sample_id.unpaired.fastq -s $chrom\_$sample_id.singleton.fastq -N


time samtools sort -n $chrom.sam | \
        samtools fastq -1 $chrom\_$sample_id.R1.fastq -2 $chrom\_$sample_id.R2.fastq -0 $chrom\_$sample_id.unpaired.fastq -s $chrom\_$sample_id.singleton.fastq -N


# .20purity.bam

#time samtools sort -n $chrom.20purity.bam | samtools fastq -1 $chrom\_$sample_id.R1.fastq -2 $chrom\_$sample_id.R2.fastq -0 $chrom\_$sample_id.unpaired.fastq -s $chrom\_$sample_id.singleton.fastq -N

#
#
#
####
# when using simulate_normal_seq
# change the commands in this way
#cd $sample_id
#pwd
#printf "samtools view -h -f 0x2 $chrom.sam >  $chrom\_paired.sam""\n"
#time samtools view -h -f 0x2 $chrom.sam >  $chrom\_$sample_id\_paired.sam
#printf "samtools fastq -N $chrom\_paired.sam -1 $sample_id/$chrom\_$sample_id.R1.fastq -2 $sample_id/$chrom\_$sample_id.R2.fastq""\n"
#time samtools view -b -h -f 0x2 $chrom.sam > $chrom\_paired.bam
#samtools reheader -c 'sed 's/SM:normal_sample/SM:SPN01_normal_sample/'' $chrom\_paired.bam
#time samtools fastq -N $chrom\_paired.sam -1 $chrom\_$sample_id.R1.fastq -2 $chrom\_$sample_id.R2.fastq
#rm $chrom\_$sample_id\_paired.sam

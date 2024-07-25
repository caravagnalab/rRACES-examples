#!/bin/bash 
#SBATCH --partition=EPYC
#SBATCH -A cdslab
#SBATCH --job-name=downsamplesam
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --output=downsampling.out
#SBATCH --error=error

# samtools view -b /orfeo/cephfs/scratch/cdslab/shared/reports_sarek/mapped/Sample_1/Sample_1.sorted.cram \
#      -T /orfeo/LTS/CDSLab/LT_storage/ref_genomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
#      -o /orfeo/cephfs/scratch/area/vgazziero/CDSlab/rRaces/downsampling/test/Sample_1.sorted.bam

# samtools index /orfeo/cephfs/scratch/area/vgazziero/CDSlab/rRaces/downsampling/test/Sample_1.sorted.bam /orfeo/cephfs/scratch/area/vgazziero/CDSlab/rRaces/downsampling/test/Sample_1.sorted.bam.bai

##################
############################## DOWNSAMPLING 

module load picard
java -jar $picard DownsampleSam \
    --INPUT ${1} \
    --OUTPUT ${2} \
    -P ${3} \
    --CREATE_INDEX \
    --REFERENCE_SEQUENCE ${4}


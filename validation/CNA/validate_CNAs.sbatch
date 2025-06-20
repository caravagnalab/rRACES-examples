#!/bin/bash
#SBATCH --job-name=validate_CNA
#SBATCH --mail-type=ALL
#SBATCH --partition=EPYC
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=validate_CNA_%j.log
#SBATCH --mem-per-cpu=1024M

# Preprocess seqz
patient=SPN03
sample=SPN03_1.1
coverage=50x
purity=0.9p
seqz=/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${patient}/sarek/${coverage}_${purity}/variant_calling/sequenza/${sample}_vs_normal_sample/${sample}_vs_normal_sample_concat.seqz.gz
cp $seqz .
gzip -d ${sample}_vs_normal_sample_concat.seqz.gz
awk '{if ( ($5 > 10) && ($4 > 3) && ($9 == "het")) {print}}' ${sample}_vs_normal_sample_concat.seqz > filtered.seqz
rm ${sample}_vs_normal_sample_concat.seqz

# Validate CNA
module load R
Rscript Validate_CNA_calls.R --sample_id ${sample} --purity ${purity} --coverage ${coverage}  


#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 50g
#SBATCH --time=2:00:00
#SBATCH --output=split.out

SPN="SPN01"
cov="100"
purity="0.6"
tool="freebayes"
outdir="/orfeo/LTS/LADE/LT_storage/lvaleriani/races/validation_data"

if [[ "$tool" == "haplotypecaller" ]]; then
    vcf="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${SPN}/sarek/${cov}x_${purity}p/variant_calling/${tool}/normal_sample/normal_sample.${tool}.filtered.vcf.gz"
elif [[ "$tool" == "freebayes" ]]; then
    vcf="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${SPN}/sarek/${cov}x_${purity}p/variant_calling/${tool}/normal_sample/normal_sample.${tool}.vcf.gz"
fi

out="${out}/$SPN/normal_sample/${tool}/vcf/"
mkdir -p $out

for i in {1..22} X Y MT
do
	bcftools view ${vcf} --regions chr${i} -o ${out}/chr${i}_normal_sample.${tool}.vcf.gz -Oz
done

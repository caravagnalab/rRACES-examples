#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 50g
#SBATCH --time=2:00:00
#SBATCH --output=split.out

# define as your
SPN="SPN03"
cov="100"
purity="0.6"

# keep as it is
outdir="/orfeo/scratch/cdslab/shared/SCOUT/${SPN}/validation/germline/vcf"
mkdir -p $outdir


for i in {1..22} X Y MT
do
  tool="haplotypecaller"
  vcf="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${SPN}/sarek/${cov}x_${purity}p/variant_calling/${tool}/normal_sample/normal_sample.${tool}.filtered.vcf.gz"
	bcftools view ${vcf} --regions chr${i} -o ${outdir}/chr${i}_normal_sample.${tool}.vcf.gz -Oz
	
	tool="freebayes"
	vcf="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${SPN}/sarek/${cov}x_${purity}p/variant_calling/${tool}/normal_sample/normal_sample.${tool}.vcf.gz"
	bcftools view ${vcf} --regions chr${i} -o ${outdir}/chr${i}_normal_sample.${tool}.vcf.gz -Oz
done

#!/bin/bash
#SBATCH --job-name=sample_processing
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --output=downsampling_%A_%a.out
#SBATCH --error=downsampling_error_%A_%a.err
#SBATCH --account=cdslab

module load samtools 

sample_list="$1"
chr="$2"
races_tmp="$3"
original_depth_tumor="$4"
original_depth_normal="$5"
purity="$6"
new_depth_tumor="$7"
frac=($(echo $purity | cut -f 2 -d "."))
new_id="$8"
new_bam="$9"
tmp_n=${races_tmp}/normal
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $sample_list)
echo "Eseguendo il job: $job"

# defining the variables
sorted_file=${races_tmp}/${sample}/splitted_sorted_bam/${chr}_${sample}.sorted.bam
new_id_out=${races_tmp}/${sample}/new_id/${chr}_${sample}.sorted.ID.bam
sam_splitted=${races_tmp}/${sample}/splitted_sorted_bam/${chr}_${sample}.bam

# sorting
if [ -f ${sorted_file} ]; then
    echo "Sorted BAM files for ${chr} and ${sample} already exist, skipping samtools sort."
else
    # Run samtools sort if the representative sorted BAM file does not exist
    samtools sort -o ${sorted_file} -T ${races_tmp}/samtools_sort_tmp ${sam_splitted}
    echo "BAM files for ${chr} and ${sample} sorted"
fi

# ID replacing
if [ -f ${new_id_out} ]; then
	echo "BAM file with new ID for ${chr} and ${sample} already exists, skipping samtools addreplacerg."
else
   # Run samtools addreplacerg if the BAM file with the new ID does not exist
   samtools addreplacerg \
    -r ID:${new_id} \
    -r SM:${new_id} \
    -r PL:ILLUMINA \
    -o ${new_id_out} \
    ${sorted_file} ## required to have the same ID once merged
  echo "id changed for ${chr} and ${sample}" 
fi

# Run downsampling
bash downsampling.sh \
    -t ${new_id_out} \
    -n ${tmp_n}/new_id/${chr}.sorted.ID.bam \
    -i $original_depth_tumor \
    -j $original_depth_normal \
    -p $purity \
    -d $new_depth_tumor \
    -o ${new_bam}/${sample}/${chr}_${sample}.purity_${frac}_coverage_${new_depth_tumor}.bam
mkdir -p ${races_tmp}/${sample}/bam2fastq
samtools sort -n ${new_bam}/${sample}/${chr}_${sample}.purity_${frac}_coverage_${new_depth_tumor}.bam -T ${races_tmp}/samtools_sort_tmp | \
	samtools fastq -1 ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.R1.fastq \
          -2 ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.R2.fastq \
          -0 ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.unpaired.fastq -s ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.singleton.fastq -N

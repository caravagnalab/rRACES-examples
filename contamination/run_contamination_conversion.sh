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

# define the folders 
# sequencing folder
races_dir="path/to/dir/where/sam/from/races/tumor/are/sequenced"
# downsampling folder
downsampling=${races_dir}/downsampling
# tmp folder
races_tmp=${downsampling}/tmp/
tmp_n=${races_tmp}normal
# normal folder 
races_dir_n="path/to/dir/where/sam/from/races/normal/are/sequenced/"

sam_files=($(ls $races_dir/*.sam | grep chr | rev | cut -f 1 -d "/"| rev |cut -f 1 -d "."))
chr=${sam_files[$SLURM_ARRAY_TASK_ID-1]} ## split by chromosomes

# get sample names
samples=$(samtools view -H ${races_dir}/${chr}.sam | grep '^@RG' | sed -n 's/.*\tSM:\([^\t]*\).*/\1/p' | sort | uniq)
# create folder for each sample in tmp 

samples_tmp_folder=$(echo $samples | sed "s@[^ ]*@${races_tmp}&@g" | sed "s@\([^ ]*\)@\1/splitted_sorted_bam@g")
mkdir -p ${samples_tmp_folder}

samples_tmp_folder_new_id=$(echo $samples | sed "s@[^ ]*@${races_tmp}&@g" | sed "s@\([^ ]*\)@\1/new_id@g")
mkdir -p ${samples_tmp_folder_new_id}

# split sam by chromosome per sample
samtools split -f ${races_tmp}/%\!/splitted_sorted_bam/%*_%\!.bam --threads 8 ${races_dir}/${chr}.sam 

# define the purity and final coverage you want to get
purity=0.20
frac=($(echo $purity | cut -f 2 -d "."))
original_depth_tumor=200
original_depth_normal=160
new_depth_tumor=100

new_bam=${downsampling}/purity_${frac}_coverage_${new_depth_tumor}
new_bam_samples=$(echo $samples | sed "s@[^ ]*@${new_bam}/&@g")
mkdir -p $new_bam_samples

# sort sam file and convert them to bam

# sort normal sample 
mkdir  ${races_tmp}/samtools_sort_tmp
samtools sort $races_dir_n/${chr}.sam -T ${races_tmp}/samtools_sort_tmp -o ${tmp_n}/splitted_sorted_bam/${chr}.sorted.bam

new_id="choose_new_id"

# modify the ID before merging
mkdir -p ${tmp_n}/new_id
samtools addreplacerg \
    -r ID:${new_id} \
    -r SM:${new_id} \
    -r PL:ILLUMINA \
    -o ${tmp_n}/new_id/${chr}.sorted.ID.bam \
    ${tmp_n}/splitted_sorted_bam/${chr}.sorted.bam ## required to have the same ID once merged

# sort tumor sample

for sample in ${samples[@]}; do
    
    # defining the variables
    sorted_file=${races_tmp}/${sample}/splitted_sorted_bam/${chr}_${sample}.sorted.bam
    new_id_out=${races_tmp}/${sample}/new_id/${chr}_${sample}.sorted.ID.bam
    sam_splitted=${races_tmp}/${sample}/splitted_sorted_bam/${chr}_${sample}.bam
    
    # sorting
    samtools sort -o ${sorted_file} -T ${races_tmp}/samtools_sort_tmp ${sam_splitted}
    
    # ID replacing
    samtools addreplacerg \
        -r ID:${new_id} \
        -r SM:${new_id} \
        -r PL:ILLUMINA \
        -o ${new_id_out} \
        ${sorted_file} ## required to have the same ID once merged
    
    /orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/contamination/downsampling.sh \
        -t ${new_id_out} \
        -n ${tmp_n}/new_id/${chr}.sorted.ID.bam \
        -i $original_depth_tumor \
        -j $original_depth_normal \
        -p $purity \
        -d $new_depth_tumor \
        -o $races_dir/${new_bam}/${sample}/${chr}_${sample}.purity_${frac}_coverage_${new_depth_tumor}.bam

    mkdir -p ${races_tmp}/${sample}/bam2fastq
	samtools fastq \
        -1 ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.R1.fastq \
        -2 ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.R2.fastq \
        -0 ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.unpaired.fastq \
        -s ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.singleton.fastq \
        -N
done 

# samtools sort $races_dir/${chr}.sam > $races_dir/${chr}.sorted.bam




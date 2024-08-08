#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=contamination_downsampling
#SBATCH -A cdslab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --output=contamination_%A_%a.out
#SBATCH --error=contamination_%A_%a.err
#SBATCH --array=12-13

module load samtools

## Set up info

# define the folders 
# sequencing folder

races_dir="$1" #"/orfeo/cephfs/scratch/area/vgazziero/CDSlab/rRaces/test/contamination_downsampling_final/simulation/sequencing_200X_basic_error_paired_350_1tumor_new_1"
races_dir_n="$2" #"/orfeo/cephfs/scratch/area/vgazziero/CDSlab/rRaces/test/contamination_downsampling_final/simulation/sequencing_160X_basic_error_paired_350_1normal_new_1_germline"


# races_dir="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/CHECK_PURITY/sequencing_100X_basic_error_paired_350_1tumor_new_1" #"path/to/dir/where/sam/from/races/tumor/are/sequenced"
# downsampling folder
downsampling=${races_dir}/downsampling
# tmp folder
races_tmp=${downsampling}/tmp/
tmp_n=${races_tmp}/normal
# normal folder 
# races_dir_n="/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/CHECK_PURITY/sequencing_80X_basic_error_paired_350_1normal_new_1" #"path/to/dir/where/sam/from/races/normal/are/sequenced/"


sam_files=($(ls $races_dir/*.sam | grep chr | rev | cut -f 1 -d "/"| rev |cut -f 1 -d "."))
chr=${sam_files[$SLURM_ARRAY_TASK_ID]} ## split by chromosomes

# get sample names
samples=$(samtools view -H ${races_dir}/${chr}.sam | grep '^@RG' | sed -n 's/.*\tSM:\([^\t]*\).*/\1/p' | sort | uniq)

# create file to iterate later 
mkdir -p ${races_tmp}${chr}
printf "%s\n" "${samples[@]}" > ${races_tmp}${chr}/"sample_names.txt"

# create folder for each sample in tmp 

samples_tmp_folder=$(echo $samples | sed "s@[^ ]*@${races_tmp}&@g" | sed "s@\([^ ]*\)@\1/splitted_sorted_bam@g")
mkdir -p ${samples_tmp_folder}

samples_tmp_folder_new_id=$(echo $samples | sed "s@[^ ]*@${races_tmp}&@g" | sed "s@\([^ ]*\)@\1/new_id@g")
mkdir -p ${samples_tmp_folder_new_id}

# split sam by chromosome per sample
# Check if a representative split BAM file already exists
#chr_21_SPN01_Sample_1.bam
if [ -f ${races_tmp}/splitted_sorted_bam/${chr}*.bam ]; then
    echo "Split BAM files for ${chr} already exist, skipping samtools split."
else
    # Run samtools split if the representative split BAM file does not exist
    samtools split -f ${races_tmp}/%\!/splitted_sorted_bam/%*_%\!.bam --threads 8 ${races_dir}/${chr}.sam
    echo "BAM files for ${chr} split"
fi

#samtools split -f ${races_tmp}/%\!/splitted_sorted_bam/%*_%\!.bam --threads 8 ${races_dir}/${chr}.sam 
#echo "splitted sam"
# define the purity and final coverage you want to get
purity="$3"
frac=($(echo $purity | cut -f 2 -d "."))
original_depth_tumor="$4"
original_depth_normal="$5"
new_depth_tumor="$6"

new_bam=${downsampling}/purity_${frac}_coverage_${new_depth_tumor}
new_bam_samples=$(echo $samples | sed "s@[^ ]*@${new_bam}/&@g")
mkdir -p $new_bam_samples

# sort sam file and convert them to bam

# sort normal sample 
mkdir  -p ${races_tmp}/samtools_sort_tmp
mkdir -p ${tmp_n}/splitted_sorted_bam
# Check if the sorted BAM file already exists
if [ -f ${tmp_n}/splitted_sorted_bam/${chr}.sorted.bam ]; then
    echo "Sorted BAM file for ${chr} already exists, skipping samtools sort."
else
    # Run samtools sort if the sorted BAM file does not exist
    samtools sort $races_dir_n/${chr}.sam -T ${races_tmp}/samtools_sort_tmp -o ${tmp_n}/splitted_sorted_bam/${chr}.sorted.bam
    echo "bam sorted for ${chr}"
fi

#samtools sort $races_dir_n/${chr}.sam -T ${races_tmp}/samtools_sort_tmp -o ${tmp_n}/splitted_sorted_bam/${chr}.sorted.bam
#echo "bam sorted"
#
new_id="$7" #"SPN01_downsampled_07P" "choose_new_id"

# modify the ID before merging
mkdir -p ${tmp_n}/new_id

# Check if the BAM file with the new ID already exists
if [ -f ${tmp_n}/new_id/${chr}.sorted.ID.bam ]; then
    echo "BAM file with new ID for ${chr} already exists, skipping samtools addreplacerg."
else
    # Run samtools addreplacerg if the BAM file with the new ID does not exist
    samtools addreplacerg \
        -r ID:${new_id} \
        -r SM:${new_id} \
        -r PL:ILLUMINA \
        -o ${tmp_n}/new_id/${chr}.sorted.ID.bam \
        ${tmp_n}/splitted_sorted_bam/${chr}.sorted.bam
    echo "id normal changed for ${chr}"
fi


#samtools addreplacerg \
#    -r ID:${new_id} \
#    -r SM:${new_id} \
#    -r PL:ILLUMINA \
#    -o ${tmp_n}/new_id/${chr}.sorted.ID.bam \
#    ${tmp_n}/splitted_sorted_bam/${chr}.sorted.bam ## required to have the same ID once merged
#echo "id normal changed"

# launch processing of the tumour bam on each sample
sample_list=${races_tmp}${chr}/sample_names.txt
n_sample_job=$(wc -l < ${sample_list})
sbatch --array=1-$(($n_sample_job)) tumour_processing.sh ${sample_list} ${chr} ${races_tmp} ${original_depth_tumor} ${original_depth_normal} ${purity} ${new_depth_tumor} ${new_id} ${new_bam}

## sort tumor sample

# for sample in ${samples[@]}; do
    
#     # defining the variables
#     sorted_file=${races_tmp}/${sample}/splitted_sorted_bam/${chr}_${sample}.sorted.bam
#     new_id_out=${races_tmp}/${sample}/new_id/${chr}_${sample}.sorted.ID.bam
#     sam_splitted=${races_tmp}/${sample}/splitted_sorted_bam/${chr}_${sample}.bam
    
#     # sorting
#     samtools sort -o ${sorted_file} -T ${races_tmp}/samtools_sort_tmp ${sam_splitted}
#     echo "tumor bam sorted"
#     # ID replacing
#     samtools addreplacerg \
#         -r ID:${new_id} \
#         -r SM:${new_id} \
#         -r PL:ILLUMINA \
#         -o ${new_id_out} \
#         ${sorted_file} ## required to have the same ID once merged
#     echo "tumor bam id changed"
#     /orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/contamination/downsampling.sh \
#         -t ${new_id_out} \
#         -n ${tmp_n}/new_id/${chr}.sorted.ID.bam \
#         -i $original_depth_tumor \
#         -j $original_depth_normal \
#         -p $purity \
#         -d $new_depth_tumor \
#         -o ${new_bam}/${sample}/${chr}_${sample}.purity_${frac}_coverage_${new_depth_tumor}.bam

#     mkdir -p ${races_tmp}/${sample}/bam2fastq
#     samtools fastq -1 ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.R1.fastq \
# 	    -2 ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.R2.fastq \
# 	    -0 ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.unpaired.fastq -s ${races_tmp}/${sample}/bam2fastq/${chr}_${sample}.singleton.fastq -N ${new_bam}/${sample}/${chr}_${sample}.purity_${frac}_coverage_${new_depth_tumor}.bam
# done

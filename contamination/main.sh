#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 --races_dir_tumor <races_dir_tumor> --races_dir_normal <races_dir_normal> --purity <purity> --depth_normal_sample <depth_normal_sample> --depth_tumor_sample <depth_tumor_sample> --new_depth_tumor_sample <new_depth_tumor_sample> --new_id <new_id>"
    exit 1
}

# Ensure the correct number of arguments are provided
if [ "$#" -ne 14 ]; then
    usage
fi

# Initialize variables
races_dir_tumor=""
races_dir_normal=""
purity=""
depth_tumor_sample=""
depth_normal_sample=""
new_depth_tumor_sample=""
new_id=""

# Parse command-line arguments
while [ "$#" -gt 0 ]; do
    case "$1" in
        --races_dir_tumor)
            races_dir_tumor="$2"
            shift 2
            ;;
        --races_dir_normal)
            races_dir_normal="$2"
            shift 2
            ;;
        --purity)
            purity="$2"
            shift 2
            ;;
        --depth_normal_sample)
            depth_normal_sample="$2"
            shift 2
            ;;
        --depth_tumor_sample)
            depth_tumor_sample="$2"
            shift 2
            ;;
        --new_depth_tumor_sample)
            new_depth_tumor_sample="$2"
            shift 2
            ;;
        --new_id)
            new_id="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            usage
            ;;
    esac
done

# Verify that all required parameters have been set
if [ -z "$races_dir_tumor" ] || [ -z "$races_dir_normal" ] || [ -z "$purity" ] || [ -z "$depth_normal_sample" ] || [ -z "$depth_tumor_sample" ] || [ -z "$new_depth_tumor_sample" ] || [ -z "$new_id" ]; then
    usage
fi

# Display the collected parameters
echo "Tumor Path: $races_dir_tumor"
echo "Normal Path: $races_dir_normal"
echo "Purity: $purity"
echo "Initial coverage of tumor samples: $depth_tumor_sample"
echo "Initial coverage of normal sample: $depth_normal_sample"
echo "New coverage of tumor samples: $new_depth_tumor_sample"

# Add the code to process the files here
# For example, we can simply list the files in the specified directories
#echo "Listing files in the tumor path: $races_dir_tumor"
#ls "$races_dir_tumor"
#
#echo "Listing files in the normal path: $races_dir_normal"
#ls "$races_dir_normal"

# Further processing can be added here
# For example, you can use the purity and depth_tumor_sample values in calculations or as parameters for other scripts/programs

# Run contamination script
sbatch run_contamination_conversion.sh "$races_dir_tumor" "$races_dir_normal" "$purity" "$depth_tumor_sample" "$depth_normal_sample" "$new_depth_tumor_sample" "$new_id" 
sample_list=${races_dir_tumor}/downsampling/tmp/chr_22/sample_names.txt
n_sample_job=$(wc -l < ${sample_list})
sbatch --array=1-$(($n_sample_job)) zcat_gzip.sh ${sample_list} "$races_dir_tumor" "$purity" "$new_depth_tumor_sample" 
exit 0


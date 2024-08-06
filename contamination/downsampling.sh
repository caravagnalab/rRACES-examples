#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -t <tumor.bam> -n <normal.bam> -i <IDP_t> -j <IDP_n> -p <purity> -d <FDP_t> -o <output.bam>"
    exit 1
}

# Parse command line arguments
while getopts "t:n:i:j:p:d:o:" opt; do
    case ${opt} in
        t) tumor_bam=${OPTARG} ;;
        n) normal_bam=${OPTARG} ;;
        i) IDP_t=${OPTARG} ;;
        j) IDP_n=${OPTARG} ;;
        p) purity=${OPTARG} ;;
        d) FDP_t=${OPTARG} ;;
        o) output_bam=${OPTARG} ;;
        *) usage ;;
    esac
done

# Check if all arguments are provided
if [ -z "$tumor_bam" ] || [ -z "$normal_bam" ] || [ -z "$IDP_t" ] || [ -z "$IDP_n" ] || [ -z "$purity" ] || [ -z "$FDP_t" ] || [ -z "$output_bam" ]; then
    usage
fi

# Ensure IDP_t, IDP_n, and FDP_t are numbers, and purity is a decimal between 0 and 1
if ! [[ $IDP_t =~ ^[0-9]+$ ]] || ! [[ $IDP_n =~ ^[0-9]+$ ]] || ! [[ $FDP_t =~ ^[0-9]+$ ]] || ! [[ $purity =~ ^0\.[0-9]+$ ]]; then
    echo "Error: IDP_t, IDP_n, and FDP_t must be integers, and purity must be a decimal between 0 and 1."
    exit 1
fi

# Calculate desired depths
D_t=$(awk "BEGIN {print $FDP_t * $purity}")
D_n=$(awk "BEGIN {print $FDP_t * (1 - $purity)}")

# Calculate downsampling fractions
f_t=$(awk "BEGIN {print $D_t / $IDP_t}")
f_n=$(awk "BEGIN {print $D_n / $IDP_n}")

# Downsample the tumor BAM file
samtools view -s $f_t -b "$tumor_bam" > "tumor_downsampled.bam"

# Downsample the normal BAM file
samtools view -s $f_n -b "$normal_bam" > "normal_downsampled.bam"

# Merge the downsampled BAM files
samtools merge -c "$output_bam" "tumor_downsampled.bam" "normal_downsampled.bam"

# Clean up temporary files
rm "tumor_downsampled.bam" "normal_downsampled.bam"

echo "Downsampling and merging completed. Output saved to $output_bam"

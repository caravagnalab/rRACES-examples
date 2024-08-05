#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input.sam> -d <IDP_t> -p <purity> -o <FDP_t>"
    exit 1
}

# Parse command line arguments
while getopts "i:d:p:o:" opt; do
    case ${opt} in
        i) input_sam=${OPTARG} ;;
        d) IDP_t=${OPTARG} ;;
        p) purity=${OPTARG} ;;
        o) FDP_t=${OPTARG} ;;
        *) usage ;;
    esac
done

# Check if all arguments are provided
if [ -z "$input_sam" ] || [ -z "$IDP_t" ] || [ -z "$purity" ] || [ -z "$FDP_t" ]; then
    usage
fi

# Calculate DP_t and f using awk
DP_t=$(awk "BEGIN {print $FDP_t * $purity}")
f=$(awk "BEGIN {print $DP_t / $IDP_t}")

# Run samtools with the calculated f
samtools view -s $f -b "$input_sam" > "${input_sam%.sam}.sorted.${purity}X.bam"

echo "Downsampling completed. Output saved to ${input_sam%.sam}.sorted.${purity}X.bam"

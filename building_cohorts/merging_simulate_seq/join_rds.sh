#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 --input_dir <input_dir> --group_size <group_size>" 
    exit 1
}

# Ensure the correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    usage
fi

# Initialize variables
input_dir=""
group_size=""


while [ "$#" -gt 0 ]; do
    case "$1" in
        --input_dir)
            input_dir="$2"
            shift 2
            ;;
        --group_size)
            group_size="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            usage
            ;;
    esac
done


#input_dir=$1
file_count=$(ls $input_dir/seq_results_SPN01_*rds | wc -l)

# 2. Define the size of each group (10 in this case)
#group_size=$2

# Initialize an empty array to store the intervals
intervals=()
chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
# 3. Loop to create the intervals
for ((i=1; i<=file_count; i+=$group_size)); do
  # Calculate the end of the interval
  end=$((i + $group_size - 1))

  # If the end exceeds the total file count, set it to the file count
  if [ "$end" -gt "$file_count" ]; then
    end=$file_count
  fi

  # Append the interval to the array
  intervals+=("$i:$end")
done

for i in ${intervals[@]}; do
  for c in ${chroms[@]}; do
    echo "sbatch -J merge_rds_$i\_$c run_merge.sh $i --lot_type single --rds_dir $input_dir --chrom $c"
    sbatch -J merge_rds_$i\_$c run_merge.sh $i single $input_dir $c
  done
done


#big_lot=$(($file_count / $group_size))
#final_interval="1:$big_lot"
#sbatch -J merge_rds_final_$final_interval run.sh $final_interval final $input_dir

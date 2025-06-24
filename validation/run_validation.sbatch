#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=validation.log
#SBATCH --job-name=validation

module load R/4.4.1

# to modify
SPN="SPN01"
DIR="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation"
# combinations you want to validate
COVERAGE=(50)
PURITY=(0.3 0.6)

# to keep as it is
QUEUE="EPYC"
ACCOUNT="cdslab"

# validate normal
echo "Submitting job for normal"
JOBID=$(sbatch --parsable \
               --job-name=validate_comb_${COV}_${PUR} \
               --partition=$QUEUE \
               --nodes=1 \
               --cpus-per-task=1 \
               --mem=24GB \
               --time=4:00:00 \
               --output=out/val_normal.out \
               --error=out/val_normal.err \
               --wrap="$DIR/validate_normal.py $SPN 50 0.3 $DIR -P $QUEUE -A $ACCOUNT")


if [ $? -eq 0 ]; then
    echo "Job submitted successfully: $JOBID"
else
    echo "Failed to submit job for COV=$COV, PUR=$PUR"
    echo "Error: $JOBID"
fi

job_ids_single_combination=()
# validate combination
for COV in "${COVERAGE[@]}"; do
  for PUR in "${PURITY[@]}"; do
    echo "Submitting job for COV=$COV, PUR=$PUR"

    JOBID=$(sbatch --parsable \
                   --job-name=validate_comb_${COV}_${PUR} \
                   --partition=$QUEUE \
                   --nodes=1 \
                   --cpus-per-task=1 \
                   --mem=24GB \
                   --time=4:00:00 \
                   --output=out/val_${COV}_${PUR}.out \
                   --error=out/val_${COV}_${PUR}.err \
                   --wrap="$DIR/validate_combination_test.py $SPN $COV $PUR $DIR -P $QUEUE -A $ACCOUNT")
    job_ids_single_combination+=("$JOBID") 
    if [ $? -eq 0 ]; then
        echo "Job submitted successfully: $JOBID"
    else
        echo "Failed to submit job for COV=$COV, PUR=$PUR"
        echo "Error: $JOBID"
    fi
  done
done

job_ids_all_combination=$(IFS=, ; echo "${job_ids_single_combination[*]}")
echo "Submitting job for $SPN"

JOBID=$(sbatch --parsable \
               --dependency=afterok:$job_ids_all_combination \
               --job-name=validate_SPN_${COV}_${PUR} \
               --partition=$QUEUE \
               --nodes=1 \
               --cpus-per-task=1 \
               --mem=24GB \
               --time=4:00:00 \
               --output=out/val_${SPN}.out \
               --error=out/val_${SPN}.err \
               # --wrap="$DIR/validate_SPN.py $SPN 50 0.3 $DIR -P $QUEUE -A $ACCOUNT")


if [ $? -eq 0 ]; then
    echo "Job submitted successfully: $JOBID"
else
    echo "Failed to submit job for SPN=$SPN"
    echo "Error: $JOBID"
fi

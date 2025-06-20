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
SPN="SPN03"
DIR="/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation"
# combinations you want to validate
COVERAGE=(50)
PURITY=(0.3 0.6)

# to keep as it is
QUEUE="EPYC"
ACCOUNT="lade"

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
                   --wrap="$DIR/validate_combination.py $SPN $COV $PUR $DIR -P $QUEUE -A $ACCOUNT")

    if [ $? -eq 0 ]; then
        echo "Job submitted successfully: $JOBID"
    else
        echo "Failed to submit job for COV=$COV, PUR=$PUR"
        echo "Error: $JOBID"
    fi
  done
done

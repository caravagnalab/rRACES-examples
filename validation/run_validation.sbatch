#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=validation.log
#SBATCH --job-name=validation

module load R/4.4.1

SPN='SPN03'
COV="100"
PUR="0.9"
DIR="/orfeo/cephfs/scratch/area/lvaleriani/ProCESS-examples"
QUEUE="GENOA"
ACCOUT="lade"

$DIR/validation/validate_combination.py $SPN $COV $PUR $DIR -P $QUEUE -A $ACCOUT

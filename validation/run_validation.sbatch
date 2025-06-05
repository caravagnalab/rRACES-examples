#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=validation.log
#SBATCH --job-name=validation

module load R/4.4.1

SPN=$1
COV=$2
PUR=$3
DIR=${PWD}

QUEUE="EPYC"
ACCOUT="cdslab"

$DIR/validation/validate_combination.py $SPN $COV $PUR $DIR -P $QUEUE -A $ACCOUT

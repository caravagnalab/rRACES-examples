#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 20gb
#SBATCH --time=96:00:00
#SBATCH --output=seq.out
#SBATCH --error=seq.err


module load singularity
user="cdslab" ## orfeo group
spn="SPN01" ## spn name
phylo="/orfeo/cephfs/fast/cdslab/ggandolfi/SPN01/phylo_forest.sff" ## path to phyloforest
out="/orfeo/cephfs/fast/cdslab/ggandolfi/SPN01/out" ## path to output directory
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/races_v3.sif" ## path to the singularity image
tmp="/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/scratch_node" ## path to a tmp directory
partition="GENOA" ## name of the partition

$path/benchmark_build_cohort.py -P $partition -A $user -s $tmp -I $image $spn $phylo $out 

#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 20gb
#SBATCH --time=96:00:00
#SBATCH --output=seq.out
#SBATCH --error=seq.err

module load singularity

# change them accordingly
partition=EPYC
user="cdslab"
spn="SPN01"

# change with your own absolute path
path="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/build_cohorts"

# keep them as they are
phylo="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/races/phylo_forest.sff"
tmp="/orfeo/cephfs/fast/cdslab/${USER}/tmp_files"
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.0.0.sif"
out="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing"
sarek_output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sarek"
tumourevo_output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/tumourevo"


$path/benchmark_build_cohort.py -P $partition -A $user -s $tmp -I $image $spn $phylo $out -C $path/orfeo.config -SD $sarek_output_dir -TD $tumourevo_output_dir

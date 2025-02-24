#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 20gb
#SBATCH --time=96:00:00
#SBATCH --output=seq.out
#SBATCH --error=seq.err

module load singularity

# change them accordingly
partition=GENOA
user="cdslab"
spn="SPN01"

# change with your own absolute path
phylo="/orfeo/cephfs/fast/cdslab/ggandolfi/SPN01/phylo_forest_NEW_SIM.sff"
tmp="/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/scripts_Alberto/scratch_node"
path="/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/scripts_Alberto"

# keep them as they are
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/races_v3.sif"
config="/orfeo/cephfs/scratch/cdslab/ggandolfi/races/sarek.config"
out="/orfeo/cephfs/fast/cdslab/shared/rRACES/${spn}/sequencing"
sarek_output_dir="/orfeo/cephfs/fast/cdslab/shared/rRACES/${spn}/sarek"

$path/build_cohort.py -P $partition -A $user -s $tmp -I $image $spn $phylo $out -C $config -SD $sarek_output_dir
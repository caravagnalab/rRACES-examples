# How to Build The Cohort from The Phylogenetic Forests

For every SPN, the cohort contains 13 sets of reads: 12 of them simulate 
the tumour sequencing with coverage in 50X, 100X, 150X, and 200X and purity
varying among 0.3, 0.6, and 0.9; the last one is the 30X simulated
sequencing of normal samples. The cells in the normal sample share the same
mutations, i.e., the germline mutations.

For efficiency purposes, not all the 12 tumour read sets were directly
produced by rRACES. Instead, for each purity, we split the generation 
of the 200X-coverage read set into 40 5X-coverage lots.
The lots are pairwise different because of their sequencing simulation
random seeds: the $i$-th lot random seed is $i$.

The 50X, 100X, 150X, 200X coverage sets were assembled by grouping 10, 20, 
30, and 40 lots, respectively.

The normal data were produced in an analogous way by splitting the computation 
of each 30X-coverage read set into lots of 5X coverage.

With the same scripts the `.rds` of the sequencing produced by `simulate_seq()` 
for each lot will be merged into a single one.

Lastly samplesheet and bash script for running sarek will be created.

<!-- The lot-based approach allowed us to distribute the read set generation
on a cluster: the more lots, the faster the generation. 
The lot generation were distributed on a computational cluster managed by SLURM.
The cluster consists of 8 nodes equipped with an AMD EPYC 7H12 64-core
processor with 512GB of RAM and 352GB of local scratch.
The number of lots, i.e., 40, is tailored on this computational facility. -->

## Generating the Tumour Cohort

This directory contains the Python scripts `build_cohort.py` that build the cohort 
of a SPN by using its phylogenetic forest file on Orfeo.

To run the scripts write a `bash` executable file like this one (we put 
`run_build_cohort.sh` with SPN01 parameters as an example):

```{sh}
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
config="/orfeo/cephfs/scratch/cdslab/ggandolfi/races/myconf.config"
out="/orfeo/cephfs/fast/cdslab/shared/rRACES/${spn}/sequencing"
sarek_output_dir="/orfeo/cephfs/fast/cdslab/shared/rRACES/${spn}/sarek"

$path/build_cohort.py -P $partition -A $user -s $tmp -I $image $spn $phylo $out -C $config -SD $sarek_output_dir
```

where you have to change: 
- `partition`: the Orfeo HPC partion where to run the scripts
- `user`: your Orfeo group 
- `spn`: the name of the SPN
- `phylo`: absolute path to the phylogenetic forest
- `tmp`: absolute path to a directory where temporaty files will be written and then deleted
- `path`: absolute path to the `build_cohort.py` script

keep as they are:
- `image`: absolute path to the singulariy image 
- `config`: absolute path to the sarek config
- `out`: the absolute path where sequencing results will be written
- `sarek_output_dir`: absolute path to the directory where sarek output will be written


Then, when you have replaced the variable, make sure that `build_cohort.py` is an 
executable file by doing:
```{sh}
chmod +x build_cohort.py
```

and finally you can sbatch the script like this:
```{sh}
sbatch run_build_cohort.sh
```

## Output files

<!-- To generate the cohort, copy the scripts into a _writable_
directory in your home on Orfeo, log in on Orfeo, and execute the following commands.

```{sh}
# clean the node local scratch directory
./clean_local_scratches.sh EPYC <USER>

# build the cohort
nohup ./build_cohort.py -P EPYC -A <USER> <SPN> <PHYLO_FOREST> <OUTPUT_DIR> &
```

`<USER>` is your SLURM user, `<SPN>` is the name of the SPN (e.g., SPN01), `<PHYLO_FOREST>` is the absolute
path of the phylogenetic forest, and `<OUTPUT_DIR>` is the absolute path of the output
directory.

For faster executions, both `<PHYLO_FOREST>` and `<OUTPUT_DIR>` should be
located in one of the `/fast` subdirectories.

> [!CAUTION]
> The content of the `/fast` subdirectories is deleted every 30 days -->

The Python script creates three directory into the `<OUTPUT_DIR>`: 
`tumour`, `normal`, and `sarek`.

The directories `tumour` and `normal` organizes the produced data by purity.
The former contains three sub-directory: `purity_0.3`, `purity_0.6`, and
`purity_0.9`. The latter has only one sub-directory, i.e., `purity_1`, meaning 
that the normal data has no contaminants.
Each purity directory stores the sequencing dataframe output (sub-directory 
`data`), the BAM file (sub-directory `BAM`), the correspoding compressed
fastq file (sub-directory `FASTQ`), and the log files of each lot.

Finally, the directory `sarek` contains X files: one for each combination of coverage 
among 50X, 100X, 150X, and 200X and purity among 0.3, 0.6, and 0.9.
There are `.csv` file with sarek datasheet for an experiment whose tumour coverage and 
purity are those declared in the file name and whose normal coverage is 30X for both 
the alignement part and the variant calling one.
In addition there will be `.sh` files that are required for running `sarek`.

```
<OUTPUT_DIR>
├── tumour                                      # tumour data
│   ├── purity_0.3                              # tumour data for purity 0.3
│   │   ├── BAM
│   │   │   ├── t00.bam                         # rRACES BAM file of lot 00
│   │   │   └── ...
│   │   ├── data
│   │   │   ├── mutations
│   │   │   │    ├── seq_results_muts_SPN01_t00.rds           # rRACES sequencing result for lot 00
│   │   │   │    └── ...
│   │   │   │    └── seq_results_muts_merged_50x_SPN01.rds    # rRACES merged sequencing results for coverage 50X     
│   │   │   │    └── ...
│   │   │   ├── parameters
│   │   │   │    ├── seq_results_params_SPN01_t00.rds         # rRACES sequencing parameters for lot 00
│   │   │   │    └── ...
│   │   │   ├── resources
│   │   │   │    ├── seq_results_resources_SPN01_t00.rds      # resources required for sequencing lot 00
│   │   │   │    └── ...
│   │   ├── FASTQ
│   │   │   ├── t00_Sample_1.R1.fastq.gz        # 5' reads of sample "Sample_1" for lot 00
│   │   │   ├── t00_Sample_1.R2.fastq.gz        # 3' reads of sample "Sample_1" for lot 00
│   │   │   ├── t00_Sample_1.singleton.fastq.gz
│   │   │   ├── t00_Sample_1.unpaired.fastq.gz
│   │   │   ├── t00_Sample_2.R1.fastq.gz        # 5' reads of sample "Sample_2" for lot 00
│   │   │   ├── t00_Sample_2.R2.fastq.gz        # 2' reads of sample "Sample_3" for lot 00
│   │   │   ├── t00_Sample_2.singleton.fastq.gz
│   │   │   ├── ...
│   │   ├── log
│   │   │   ├── lot_t00.log                     # lot 00 log 
│   │   │   └── ...
│   │   ├── t00_BAM.done                        # lot 00 BAM was produced
│   │   ├── t00_final.done                      # lot 00 fastq files were produced
│   │   └── ...
│   ├── purity_0.6                              # tumour data for purity 0.6
│   │   └── ...
│   └── purity_0.9                              # tumour data for purity 0.9
│       └── ...
├── normal                                      # normal data
│   └── purity_1                                # the normal purity is 1
│       └── ...
└── sarek
    ├── sarek_normal.csv                        # sarek datasheet for normal sample
    ├── sarek_50x_0.3p.csv                      # sarek datasheet for tumour coverage 50X and purity 0.3 
    ├── ...                     
    ├──sarek_mapping_normal.sh                 # .sh file for running sarek mapping for normal sample
    ├──sarek_mapping_50x_0.3p.sh               # .sh file for running sarek mapping for tumour coverage 50X and purity 0.3 
    ├── ...               
    ├──sarek_variant_calling_50x_0.3p.csv      # sarek datasheet for variant calling for tumour coverage 50X and purity 0.3 
    ├── ...             
    ├──sarek_variant_calling_50x_0.3p.sh       # .sh file for running sarek variant calling for tumour coverage 50X and purity 0.3 
    ├── ...             
```

## Run sarek
### Mapping
As described previously in the `sarek` directory will be stored both `.csv` and `.sh` files for running sarek steps. 
The first files that will be generated are `sarek_normal.csv` and `sarek_mapping_normal.sh`. Once these files will be written you can start launching them in this way:
```
sbatch sarek_mapping_normal.sh
```
and in this way the sarek mapping of the normal sample will start.

Then once the first purity combination (0.3) will ends, they will be written also `sarek_50x_0.3p.csv`, `sarek_100x_0.3p.csv`, `sarek_150x_0.3p.csv`, `sarek_200x_0.3p.csv`  and `sarek_mapping_50x_0.3p.sh`, `sarek_mapping_100x_0.3p.sh`,... files.

Then you can start to send the following jobs:
```
sbatch sarek_mapping_50x_0.3p.sh
sbatch sarek_mapping_100x_0.3p.sh
sbatch sarek_mapping_150x_0.3p.sh
sbatch sarek_mapping_200x_0.3p.sh
```
Pay attention that at the same time also `sarek_variant_calling_50x_0.3p.csv` and `sarek_variant_calling_50x_0.3p.csv` will be generated but at that moment they have not be runned.

Then after 0.3, will be runned 0.6 and finally 0.9 purity. After the ends of each of those you can run as before the corresponding mapping scripts.

### Variant Calling
Finally when the mapping, of a given purity, will be over you can start to run the variant calling. For example, when the mapping of 50X_0.3p will ends you can launch the following script:
```
sbatch sarek_variant_calling_50x_0.3p.sh 
``` 
and so on for each combination of coverage and purity!
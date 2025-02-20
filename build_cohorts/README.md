# How to Build The Cohort from The Phylogenetic Forests

For every SPN, the cohort contains 13 sets of reads: 12 of them simulate 
the tumour sequencing with coverage in 50X, 100X, 150X, and 200X and purity
varying among 0.3, 0.6, and 0.9; the last one is the 50X simulated
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
of each 50X-coverage read set into 40 lots.

The lot-based approach allowed us to distribute the read set generation
on a cluster: the more lots, the faster the generation. 
The lot generation were distributed on a computational cluster managed by SLURM.
The cluster consists of 8 nodes equipped with an AMD EPYC 7H12 64-core
processor with 512GB of RAM and 352GB of local scratch.
The number of lots, i.e., 40, is tailored on this computational facility.

## Generating the Tumour Cohort

This directory contains 2 Python scripts that build the cohort of a SPN by
using its phylogenetic forest file on Orfeo.

To generate the cohort, copy the scripts into a _writable_
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
> The content of the `/fast` subdirectories is deleted every 30 days

The Python script creates three directory into the `<OUTPUT_DIR>`: 
`tumour`, `normal`, and `sarek`.

The directories `tumour` and `normal` organizes the produced data by purity.
The former contains three sub-directory: `purity_0.3`, `purity_0.6`, and
`purity_0.9`. The latter has only one sub-directory, i.e., `purity_1`, meaning that the normal data has no contaminants.
Each purity directory stores the sequencing dataframe output (sub-directory 
`data`), the BAM file (sub-directory `BAM`), the correspoding compressed
fastq file (sub-directory `FASTQ`), and the log files of each lot.

Finally, the directory `sarek` contains 12 CSV files: one for each combination 
of coverage among 50X, 100X, 150X, and 200X and purity among 0.3, 0.6, and 0.9.
Each CSV file is a sarek datasheet for an experiment whose tumour coverage and 
purity are those declared in the file name and whose normal coverage is 50X.

```
<OUTPUT_DIR>
├── tumour                                      # tumour data
│   ├── purity_0.3                              # tumour data for purity 0.3
│   │   ├── BAM
│   │   │   ├── t00.bam                         # rRACES BAM file of lot 00
│   │   │   └── ...
│   │   ├── data
│   │   │   ├── seq_results_SPN01_t00.rds       # rRACES sequencing result for lot 00
│   │   │   └── ...
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
    ├── sarek_50x_0.3p.csv                      # sarek datasheet for tumour coverage 50X
    ├── ...                                     # and purity 0.3 with normal 50X
    ├── sarek_150x_0.6p.csv                     # sarek datasheet for tumour coverage 150X
    └── ...                                     # and purity 0.6 with normal 50X 
```
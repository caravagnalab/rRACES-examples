# How to Build The Cohort from The Phylogenetic Forests

For every SPN, the cohort contains 13 sets of reads: 12 of them are produced by
simulating the tumour sequencing with coverage in 50X, 100X, 150X, and 200X and
purity varying among 0.3, 0.6, and 0.9; the last one is the 30X simulated
sequencing of normal samples. The cells in the normal sample share the same
mutations, i.e., the germline mutations.

In order to have access to the simulated mutations, all 12 coverage-purity combinations are 
generated directly by using rRACES `simulate_seq` parameters of coverage and purity. To parallelize
the process the `deploy_slurm_rRACES_seq_purity.py` must be run in the follwing way:

```bash
nohup ./deploy_slurm_rRACES_seq_purity.py <SPN_ID> <path_to_phylo_forest> <path_to_output_dir> \
    -p <partition> -A <user_group> -c <coverage> -pu <purity> -l <n_of_lots> -j <n_of_parallel_jobs> & 
```

The output folder will be structured in this way:

```bash
├── BAM
├── data
├── log
├── t0.done
├── t1.done
├── t2.done
├── t3.done
├── t4.done
├── t5.done
├── t6.done
├── t7.done
├── t8.done
└── t9.done
```

In the `BAM` folder you will find the sequenced bam files while in the `data` folder you will find all the `rds` files
containing the `simulate_seq` object. If the coverage is set to 200 and the number of lots is set to 20, you will find 
20 bam files and 20 rds files representing a 20X sequencing depth.

## Merging of single lot `simulate_seq` objects

In order to have the final `simulate_seq` object recreating the initial condition of coverage and purity, we need to
merge the single `rds` files. To do so, we will split the total number of `rds` into small chunks to parallelize the merging. It is
important that the merging script takes in input no more than 10 `rds` files. First we have to merge the single `rds` into small chunks:

```bash
bash join_rds.sh --input_dir <path_to_output_dir/data> --group_size 10
```
This script will parallelize in terms if chunks and chromosomes. Once the process is finished, then run the final merge:

```bash
sbatch -J merge_rds_final run_merge.sh 1:5 final <path_to_output_dir/data>
```

The first argument is the final interval that must be merged.

### Example

Sequence at 200X in 20 lots (each lot will have coverage equals to 10X).

```bash
nohup ./deploy_slurm_rRACES_seq_purity.py <SPN_ID> <path_to_phylo_forest> <path_to_output_dir> \
    -p <partition> -A <user_group> -c 200 -pu <purity> -l 20 -j 20 &
```
The run the first merging by chunks of size 5:
```bash
bash join_rds.sh --input_dir <path_to_output_dir/data> --group_size 5
```
In this way you will end up with 4 chunks named `seq_results_<chrom>_lot_1:5.rds`, `seq_results_<chrom>_lot_6:10.rds`,
`seq_results_<chrom>_lot_11:15.rds` and `seq_results_<chrom>_lot_16:20.rds`. Now you can run the final merge:

```bash
sbatch -J merge_rds_final run_merge.sh 1:4 final <path_to_output_dir/data>
```

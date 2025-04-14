> [!WARNING]
> Outdated information. Please refer to [this directory](https://github.com/caravagnalab/ProCESS-examples/tree/main/building_cohorts) to produce the cohort.

# Downsampling and data preparation

## 1. ProCESS sequencing
In order to have different combination of purity and coverage starting from the same simulation, we need to sequence two times:

1. Run `simulate_seq` without the normal sample at a maximum coverage;
2. Run `simulate_normal_seq` at a maximum coverage required for downsampling (e.g. if tumor depth is set to 200X, the maximum coverage for the normal is 160X).

In particular, the command to sequence a tumor sample is:

```r
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose

seq_results <- parallel::mclapply(chromosomes, function(c) {
	simulate_seq(phylo_forest, chromosomes = c, coverage = 100,write_SAM = TRUE,
		     sequencer = basic_seq, insert_size_mean = 350, insert_size_stddev = 10,
		     output_dir = "/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/CHECK_PURITY/sequencing_100X_basic_error_paired_350_1tumor_new_1",
	update_SAM =TRUE, with_normal_sample =FALSE, template_name_prefix = "td")
}, mc.cores = 8)
```

and the command to sequence a normal sample is:

```r
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose
seq_results <- parallel::mclapply(chromosomes, function(c) {
    simulate_normal_seq(phylo_forest, chromosomes = c, coverage = 80,write_SAM = TRUE, read_size =150,
          sequencer = basic_seq, insert_size_mean = 350, insert_size_stddev = 10,
		  output_dir = "/orfeo/LTS/CDSLab/LT_storage/ggandolfi/races_simulations/CHECK_PURITY/sequencing_80X_basic_error_paired_350_1normal_new_1_germline",
    update_SAM =TRUE, with_preneoplastic = TRUE, template_name_prefix = "nd")
}, mc.cores = 8)
```

Once the sequencing is finished, each folder will contain:

- A `${chrom}.sam` for each chromosome. If the SPN is multisample, each `.sam` will contains all the sample reads merged togheter. So a furhter step of splitting will be necessary to have sample specific `.bam` files.
- A `${sample}/${crhom}.dat` file containg coverage information (not useful).

Everythin is now ready to perform downsampling and conversion to `fastq` files.

## 2. Contamination

Starting from an initial tumor sample at high coverage and 100% purity, we can generate different combination of coverage and purity, for example:

- change coverage, keep purity fixed;
- change purity, keep coverage fixed;
- change both purity and coverage.

In order to do so, you simply have to run the `run_contamination_conversion.sh` script which requires in input:

- the path to the `simulate_seq` output;
- the path to the `simulate_normal_seq` output;
- the original sequencing depth of tumor and normal samples;
- the expected new coverage and purity of tumor sample;
- the name of the new sample.

This script will perform:

1. Contamination of the tumor sample with their corresponding normal resulting in a new `${chrom}.bam` file. The fraction $f_{tumor}$ to downsample the tumor `bam` is derived from the intial coverage $ODP_{tumor}$, the expected new coverage $FDP_{tumor}$ and the expected new purity $\pi$:  

    $$DP_{tumor} =FDP_{tumor} × \pi$$ \
    $$f_{tumor} =\frac{DP_{tumor}}{ODP_{tumor}}$$

    While the fraction $f_{normal}$ to downsample the normal `bam`:

    $$DP_{normal} =FDP_{tumor} × (1-\pi)$$ \
    $$f_{normal} =\frac{DP_{normal}}{ODP_{normal}}$$

2. Conversion of the `bam` file to paired-end `fastq` files, required to run nf-core/sarek.

The output folder of the contamination script will be strucutred in this way:

``` bash
downsampling/
├── purity_${purity}_coverage_${coverage}
│   └── ${sample}
│       └── ${chr}_${sample}.purity_${purity}_coverage_${coverage}.bam
└── tmp
    ├── normal
    │   ├── new_id
    │   │   └── ${chr}.sorted.ID.bam
    │   └── splitted_sorted_bam
    │       └── ${chr}.sorted.bam
    ├── samtools_sort_tmp
    └── ${sample}
        ├── bam2fastq
        │   ├── ${chr}_${sample}.R1.fastq
        │   ├── ${chr}_${sample}.R2.fastq
        │   ├── ${chr}_${sample}.singleton.fastq
        │   └── ${chr}_${sample}.unpaired.fastq
        ├── new_id
        │   └── ${chr}_${sample}.sorted.ID.bam
        └── splitted_sorted_bam
            ├── ${chr}_${sample}.bam
            └── ${chr}_${sample}.sorted.bam
```

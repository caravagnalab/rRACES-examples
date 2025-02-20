# SPN01

## Tissue dynamics and Sampling

In the `0_sim_tissue.R` script the tumour dynamics for SPN01 is generated. In order to generate the final report, some mandatory steps are required:

1. Store tissue plots and state plots in different moment of the simulation that you think are important to understand the dynamics;
2. Save an `rds` containing the muller plot, since it will be used for plotting signature exposure in time.

## Setting up the `MutationEngine`

In `1_mut_engine.R` script we can add information about driver mutations and CNA, passanger mutation rates, germline sample and signature exposure to the previously saved `sample_forest.sff`. This script will also save the second part of the report in the `page2.rds` object.

## Simulate sequencing

In `2_sequencing.R` script sequencing is simulated starting from the previously saved `phylo_forest.sff`. In this script you have to set up:

1. Illumina sequencing error rate;
2. Coverage and purity;
3. General paramters for saving `bam` files. At this level, we are not interest in saving the real `bam` since a more efficient procedure wìll be used.

This script will output a dataframe containg information about the sequenced mutations. 

## Plotting mutations and CNA



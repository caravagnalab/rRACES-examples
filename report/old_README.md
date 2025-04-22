# How to build SPN report

For every SPN, the following data are necesary to generate the final report and must be saved in the shared folder at `/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN{ID}`:

1. **simultaion snapshot**: each `0_sim_tissue.R` script need to save the snapshot of the simulation, in order to be recovered by the report script. You simiply neeed to specify `sim <- SpatialSimulation("SPN01", save_snapshots=TRUE)`.
2. **samples forest**: each `0_sim_tissue.R` script need to save the samples forest of containing the phylogeny of the sampled cells at `/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN{ID}/sample_forest.sff`;
3. **phylogenetic forest**: each `1_mut_engine.R` script need to save the phylogenetic  forest of containing the phylogeny of the sampled cells and the genomic mutations occurring in the cells at `/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN{ID}/phylo_forest.sff`;
4. **sequencing results**: each SPN should follow the build_cohort documentation to generated all the 12 coverage purity combinations (+ the normal sample). For more details see ([build cohort section](https://github.com/caravagnalab/ProCESS-examples/blob/main/build_cohorts/README.md)) Please provide the path to one of the combinations and the normal sample sequencing result. The specific coverage and purity combination selected to be showed in the report must be specified in the param.yml file.
5. **segmentation tables**: each `1_mut_engine.R` script need to save the output of `phylo_forest$get_bulk_allelic_fragmentation()` for each sample as an `.rds` file in order to have access to CNA informations.

For each SPN, you need to fill a `params_SPN{ID}.yml` where you have to specify few information. The `params_SPN{ID}.yml` file is formatted as follows:

```yaml
spn: "SPN01"
sequencing:
  coverage: 100
  purity: 0.9
files:
  sim: "/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN01/SPN01"
  sample_forest: "/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN01/sample_forest.sff"
  phylo_forest: "/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN01/phylo_forest.sff"
  seq_res: "/path/to/sequencing_result_100X_0.9p.rds"
cna_dir: "/path/to/cna/dir"
```

Once you fill in the `params_SPN{ID}.yml` open an Rstudio session and run:

```r
rmarkdown::render("Report.Rmd", params = yaml::read_yaml("params_SPN{ID}.yml"),output_dir="/orfeo/cephfs/scratch/cdslab/shared/races/report")
```

The final report will be generated in the shared scratch folder and will inherit the name of the SPN, as well as the selected coverage and purity values.

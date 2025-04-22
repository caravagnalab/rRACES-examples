# How to build SPN report

## Before running the report
Before rendering the report you should run a script for creating some dataframe that are required later.

Open `run_map.sh` and modify:
```
# path to tmp_dir
dir="/orfeo/fast/area/lvaleriani/tmp_nf"
# path to ProCESS-examples repo, path to this repo
base="/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/report"
# name of your spn
spn=SPN03
```

and then 
```
sbatch run_map.sh
```


## Running the report
For every SPN, the following data are necessary to generate the final report and must be saved in the shared folder at `/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN{ID}`:

1. **simultaion snapshot** folder at `/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN{ID}/process/sample_forest.sff`
2. **samples forest** at `/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN{ID}/process/SPN{ID}`;
3. **phylogenetic forest** at `/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN{ID}/process/phylo_forest.sff`
4. **sequencing results** at `/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN{ID}/sequencing/`. For more details see ([build cohort section](https://github.com/caravagnalab/ProCESS-examples/blob/main/build_cohorts/README.md)).
5. **cna_data** at `/orfeo/cephfs/scratch/cdslab/shared/races/SCOUT/SPN{ID}/process/cna_data`

Provvisory:
Once you have all your data you should run the `render_report.R` script in Rstudio_server with the `orfeo_rstudio_v1.sif` singularity image.

You need to change only 2 variables in the script:
- `workdir  <- "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/report"` where to put the path to this directory
- `spn <- "SPN03"` where you should put your SPN name


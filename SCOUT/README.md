# SPN simulations
To ensure uniformity in the creation of each SPN, certain requirements must be met.  

## Saving Snapshot  
Remember to save the snapshot of the simulations by running:  

```{r}
sim <- SpatialSimulation(name = 'SPN01', seed = 12345, save_snapshot=TRUE)
```

## Sample naming requirments
To standardize sample names, the following naming convention must be used:
```
SPN{X}_{T}{N}
```
where:
- `X` is the SPN number (eg, 01, 02, ...)
- `T` is the sampling time (1 = first time point, 2 = second time point, etc.)
- `N` is the sample number for a given time point T (1 = first sample at T time point, 2 = second sample at T time point)

If you are working on SPN01 and are sampling your third sample at the second time point, you must name it: `SPN01_2_3`


## Required script and files
For each SPN, the following files must be present:
- `simulate_tissue.R`
- `simulate_mutation.R`
- `params.yml` (this file is described in details in `report/README.md`)

These files must be pushed to the rRACES-examples GitHub repository in the corresponding folder SCOUT/SPNX.

## Required output files
The following output files must be saved in the directory:
`/orfeo/cephfs/scratch/cdslab/shared/rRACES/SCOUT/SPNX/races` (replace X with the number of your SPN):

- `samples_forest.sff`
- `phylogenetic_forest.sff`
- `<sample>_cna_data.rds`
- snapshot folder

The `<sample>_cna.rds` has to be obtained by running the following command for each sample:
```{r}
phylo_forest$get_bulk_allelic_fragmentation()
```


## How to run Rscript with singularity image
To run the R scripts, a Singularity image with the latest version of rRACES must be used. You can do this with the following script:

```{sh}
#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 100gb
#SBATCH --time=1:00:00
#SBATCH --output=rRACES.out
#SBATCH --error=rRACES.err

module load singularity

image='/orfeo/cephfs/scratch/cdslab/shared/SCOUT/races_v4.sif'

singularity exec --bind /orfeo:/orfeo --no-home $image Rscript simulate_tissue.R

singularity exec --bind /orfeo:/orfeo --no-home $image Rscript simulate_mutation.R
```

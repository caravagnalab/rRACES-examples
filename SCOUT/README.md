# SPN simulations
To ensure uniformity in the creation of each SPN, certain requirements must be met.  

## Simulation parameters
### Number of cells to sample 
Remember to sample around **2k cells** per sample.

### Mutations rates
Set mutation rates in MutationEngine to:
- SNV: `1e-8` (except for the SPN02 that is hypermutant)
- ID: `1e-9`
- CNA: depends on tumour type (eg. SPN01 simulation have a lot of CN events and 1e-10 rate, low number of CN with 1e-12 rate)

### Set tumour type in MutationEngine
```
# check available tumour type
get_available_tumours_in(setup_code = "GRCh38") %>% head()


m_engine <- MutationEngine(setup_code = "GRCh38",
                           tumour_type = "COAD",
                           tumour_study = "US")
```

### Set gender of simulation
```
#return available germline subject with gender
m_engine$get_germline_subjects()

#select the one that you want based on gender that you need
m_engine$set_germline_subject("")
```

Male samples:
- `SPN01`
- `SPN03`
- `SPN04`
- `SPN06`

Female samples:
- `SPN02`
- `SPN05`
- `SPN07`

### Set pre-neoplastic mutations
```
phylo_forest <- m_engine$place_mutations(samples_forest, num_of_preneoplatic_SNVs = 800, num_of_preneoplatic_indels = 200)
```


## Saving snapshot and set seed
Remember to save the snapshot of the simulations and set the seed by running:  

```{r}
set.seed(12345)
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

These files must be pushed to the ProCESS-examples GitHub repository in the corresponding folder SCOUT/SPNX.

## Required output files
The following output files must be saved in the directory:
`/orfeo/cephfs/scratch/cdslab/shared/ProCESS/SCOUT/SPNX/races` (replace X with the number of your SPN):

- `samples_forest.sff`
- `phylogenetic_forest.sff`
- `<sample>_cna.rds`
- snapshot folder

The `<sample>_cna.rds` has to be obtained by running the following command for each sample:
```{r}
phylo_forest$get_bulk_allelic_fragmentation()
```

Eg.
```
sample_names <- phylo_forest$get_samples_info()[["name"]]
lapply(sample_names,function(s){
    cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
    saveRDS(file=paste0(outdir,"cna_data/",s,"_cna.rds"),object=cna)
})
```

## How to run Rscript with singularity image
To run the R scripts, a Singularity image with the latest version of ProCESS must be used. You can do this with the following script:

```{sh}
#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 50gb
#SBATCH --time=1:00:00
#SBATCH --output=ProCESS.out
#SBATCH --error=ProCESS.err

module load singularity
image='/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.0.0.sif'

singularity exec --bind /orfeo:/orfeo --no-home $image Rscript simulate_tissue.R
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript simulate_mutation.R
```

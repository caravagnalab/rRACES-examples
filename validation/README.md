# How to run mutations validation
## Required Rpackages
For running the validation the following R packages must be installed for R version `R/4.4.1`

- `optparse`
- `tidyverse`
- `vcfR`
- `caret`
- `patchwork`
- `RColorBrewer`
- `ggExtra`
- `ggVennDiagram`
- `ggalluvial`
- `ggplotify`
- `ggplot2`
- `ProCESS`
- `CNAqc`
- `data.table`


## How to run
In order to run the validation of CNA, Somatic and Germline mutations you must run the `run_validation.sbatch` script in this way:

```
sbatch run_validation.sbatch <SPN> <COV> <PUR>
```

substituting SPN, COV and PUR with the given combination you want to validate.

Eg.
```
sbatch run_validation.sbatch SPN03 50 0.3
```






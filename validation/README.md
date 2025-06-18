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
In order to run the validation of CNA, Somatic and Germline mutations for the first time, you must open the `run_validation.sbatch` and modify the following parts:

```
# to modify
SPN="SPN03"
DIR="/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/validation"

# combinations you want to validate
COVERAGE=(50)
PURITY=(0.3 0.6)
```
substituting SPN, COV and PUR with the given combination you want to validate.


 <!-- In case you would like to run only some steps of the validation, you must specify which step to skip. Possible values are: `cna` (to skip rerunning CNA report), `somatic` (to skip rerunning somatic validation) and `germline` (to skip rerunning germline validation). Multiple values can be specified via commas. For example, if you want to run validation report only for CNAs, you simply have to run:

Eg.
```
sbatch run_validation.sbatch SPN03 50 0.3 somatic,germline
```
-->






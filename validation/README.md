# How to run mutations validation
In order to run the validation of CNA, Somatic and Germline mutations you must modifiy the `run_validation.sbatch` script in this way:
- `SPN="SPN03"`: the id of your SPN
- `COV="100"`: coverage for which you want to run the validation (sarek variant calling result)
- `PUR="0.9"`: coverage for which you want to run the validation (sarek variant calling result)
- `DIR="/orfeo/cephfs/scratch/area/lvaleriani/ProCESS-examples"`: path to `ProCESS-examples` directory

and then:
```
sbatch run_validation.sbatch
```




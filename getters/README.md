# How getter works

- [ProCESS](#process_getters)
- [sarek](#sarek_getters)
- [tumourevo](#tumourevo_getters)

> [!IMPORTANT]
> If you need some files than are not returned please let us know (open a pull request or write on slack)!


> [!IMPORTANT]
> Please do not modify the functions but let us know if there are problems or improvement to be done!


## process_getters
`ProCESS` getters return the path to the specific file/infomation.
```
> get_sample_names(spn = 'SPN04')
[1] "SPN04_1.1" "SPN04_1.2"

> get_process_cna(spn = 'SPN04')
$SPN04_1.1
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN04/process/cna_data/SPN04_1.1_cna.rds"

$SPN04_1.2
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN04/process/cna_data/SPN04_1.2_cna.rds"

> get_process_cna(spn = 'SPN04', sample = 'SPN04_1.1')
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN04/process/cna_data/SPN04_1.1_cna.rds"

> get_process_gender(spn = 'SPN02')
[1] "XX"

> get_phylo_forest(spn = 'SPN01')
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/process/phylo_forest.sff"

> get_sample_forest(spn = 'SPN01')
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/process/sample_forest.sff"

> get_mutations(spn = 'SPN01', type = 'normal')
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sequencing/normal/purity_1/data/mutations/seq_results_muts_merged_coverage_30x.rds"

> get_mutations(spn = 'SPN01', type = 'tumour', coverage =100, purity = 0.3)
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sequencing/tumour/purity_0.3/data/mutations/seq_results_muts_merged_coverage_100x.rds"
```

## sarek_getters
`sarek` getters return a list of files corresponding to the results of a combination of: 
- SPN 
- coverage
- purity
- caller
- type - when available
- sampleID - when available

### get_sarek_cna_file

```
> get_sarek_cna_file(spn = 'SPN01', 
+                    coverage = 50, 
+                    purity = 0.3, 
+                    caller = 'ascat', 
+                    type = "tumour", 
+                    sampleID = "SPN01_1.1")
$cnvs
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/ascat/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.cnvs.txt"

$purityploidy
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/ascat/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.purityploidy.txt"

$segments
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/ascat/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.segments.txt"

$tumourBAF
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/ascat/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.tumour_tumourBAF.txt"

$tumourLogR
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/ascat/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.tumour_tumourLogR.txt"

> get_sarek_cna_file(spn = 'SPN01', 
+                    coverage = 50, 
+                    purity = 0.3, 
+                    caller = 'cnvkit', 
+                    type = "tumour", 
+                    sampleID = "SPN01_1.1")
$cnr
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/cnvkit/SPN01_1.1_vs_normal_sample/SPN01_1.1.cnr"

$somatic.call
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/cnvkit/SPN01_1.1_vs_normal_sample/SPN01_1.1.somatic.call.cns"

> get_sarek_cna_file(spn = 'SPN03', 
+                    coverage = 50, 
+                    purity = 0.3, 
+                    caller = 'sequenza', 
+                    type = "tumour", 
+                    sampleID = "SPN03_1.1")
$confints_CP
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN03/sarek/50x_0.3p/variant_calling/sequenza/SPN03_1.1_vs_normal_sample/SPN03_1.1_vs_normal_sample_confints_CP.txt"

$mutations
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN03/sarek/50x_0.3p/variant_calling/sequenza/SPN03_1.1_vs_normal_sample/SPN03_1.1_vs_normal_sample_mutations.txt"

$segments
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN03/sarek/50x_0.3p/variant_calling/sequenza/SPN03_1.1_vs_normal_sample/SPN03_1.1_vs_normal_sample_segments.txt"

```

### get_sarek_vcf_file
```
> get_sarek_vcf_file(spn = 'SPN01', 
+                    coverage = 50, 
+                    purity = 0.3,
+                    caller = 'mutect2', 
+                    type = "tumour")
$vcf
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/mutect2/SPN01/SPN01.mutect2.filtered.vcf.gz"

$tbi
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/mutect2/SPN01/SPN01.mutect2.filtered.vcf.gz.tbi"

> get_sarek_vcf_file(spn = 'SPN01', 
+                    coverage = 50, 
+                    purity = 0.3,
+                    caller = 'freebayes', 
+                    type = "normal", 
+                    sampleID = 'SPN01_1.1')
$vcf
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/freebayes/normal_sample/normal_sample.freebayes.vcf.gz"

$tbi
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/freebayes/normal_sample/normal_sample.freebayes.vcf.gz.tbi"

> get_sarek_vcf_file(spn = 'SPN01', 
+                    coverage = 50, 
+                    purity = 0.3,
+                    caller = 'freebayes', 
+                    type = "tumour", 
+                    sampleID = 'SPN01_1.1')
$vcf
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/freebayes/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.freebayes.vcf.gz"

$tbi
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/freebayes/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.freebayes.vcf.gz.tbi"

> get_sarek_vcf_file(spn = 'SPN01', 
+                    coverage = 50, 
+                    purity = 0.3,
+                    caller = "haplotypecaller", 
+                    type = "normal", 
+                    sampleID = 'SPN01_1.1')
$vcf
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/haplotypecaller/normal_sample/normal_sample.haplotypecaller.filtered.vcf.gz"

$tbi
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/haplotypecaller/normal_sample/normal_sample.haplotypecaller.filtered.vcf.gz.tbi"

> get_sarek_vcf_file(spn = 'SPN01', 
+                    coverage = 50, 
+                    purity = 0.3,
+                    caller = 'strelka', 
+                    type = "tumour", 
+                    sampleID = 'SPN01_1.1')
$indels_vcf
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/strelka/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.strelka.somatic_indels.vcf.gz"

$indels_tbi
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/strelka/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.strelka.somatic_indels.vcf.gz.tbi"

$snvs_vcf
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/strelka/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.strelka.somatic_snvs.vcf.gz"

$snvs_tbi
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT//SPN01/sarek/50x_0.3p/variant_calling/strelka/SPN01_1.1_vs_normal_sample/SPN01_1.1_vs_normal_sample.strelka.somatic_snvs.vcf.gz.tbi"

```

## tumourevo_getters
`tumourevo` getters return a list of files corresponding to the results of a combination of: 
- spn 
- coverage
- purity
- vcf_caller
- cna_caller
- sample


```
## Driver ####
> get_tumourevo_driver(spn = 'SPN03', 
+                      coverage = 50, 
+                      purity = 0.3, 
+                      vcf_caller = 'mutect2', 
+                      cna_caller = 'ascat',
+                      sample = 'SPN03_1.1')
[1] "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN03/tumourevo/50x_0.3p_mutect2_ascat/driver_annotation/annotate_driver/SCOUT/SPN03/SPN03_SPN03_1.1/SCOUT_SPN03_SPN03_SPN03_1.1_driver.rds"


## Subclonal ####
## tool: mobster, pyclonevi, viber, ctree
get_tumourevo_subclonal(spn = 'SPN03', 
                     coverage = 50, 
                     purity = 0.3, 
                     vcf_caller = 'mutect2', 
                     cna_caller = 'ascat',
                     sample = 'SPN03_1.1', 
                     tool = 'mobster')

get_tumourevo_subclonal(spn = 'SPN03', 
                        coverage = 50, 
                        purity = 0.3, 
                        vcf_caller = 'mutect2', 
                        cna_caller = 'ascat',
                        sample = 'SPN03_1.1', 
                        tool = 'pyclonevi')

get_tumourevo_subclonal(spn = 'SPN03', 
                        coverage = 50, 
                        purity = 0.3, 
                        vcf_caller = 'mutect2', 
                        cna_caller = 'ascat',
                        sample = 'SPN03_1.1', 
                        tool = 'viber')

# sample specific ctree result
get_tumourevo_subclonal(spn = 'SPN03', 
                        coverage = 50, 
                        purity = 0.3, 
                        vcf_caller = 'mutect2', 
                        cna_caller = 'ascat',
                        sample = 'SPN03_1.1', 
                        tool = 'ctree')

# patient specific ctree result
get_tumourevo_subclonal(spn = 'SPN03', 
                        coverage = 50, 
                        purity = 0.3, 
                        vcf_caller = 'mutect2', 
                        cna_caller = 'ascat',
                        sample = 'SPN03', 
                        tool = 'ctree')

## QC ####
## tool: CNAqc, tinc, join_CNAqc
get_tumourevo_qc(spn = 'SPN03', 
                 coverage = 50, 
                 purity = 0.3, 
                 vcf_caller = 'mutect2', 
                 cna_caller = 'ascat',
                 sample = 'SPN03_1.1',
                 tool = 'CNAqc')

get_tumourevo_qc(spn = 'SPN03', 
                 coverage = 50, 
                 purity = 0.3, 
                 vcf_caller = 'mutect2', 
                 cna_caller = 'ascat',
                 sample = 'SPN03_1.1',
                 tool = 'tinc')

get_tumourevo_qc(spn = 'SPN03', 
                 coverage = 50, 
                 purity = 0.3, 
                 vcf_caller = 'mutect2', 
                 cna_caller = 'ascat',
                 sample = 'SPN03_1.1',
                 tool = 'join_CNAqc')

## Signature ##
## tool: SparseSignatures, SigProfiler

get_tumourevo_signatures(spn = 'SPN03', 
                 coverage = 50, 
                 purity = 0.3, 
                 vcf_caller = 'mutect2', 
                 cna_caller = 'ascat',
                 tool = 'SparseSignatures')


get_tumourevo_signatures(spn = 'SPN03', 
                         coverage = 50, 
                         purity = 0.3, 
                         vcf_caller = 'mutect2', 
                         cna_caller = 'ascat',
                         tool = 'SigProfiler',
                         context = 'SBS96') ## specify also the context

```

# How to Build The Cohort from The Phylogenetic Forests

For every SPN, the cohort contains 13 sets of reads: 12 of them are produced by
simulating the tumour sequencing with coverage in 50X, 100X, 150X, and 200X and
purity varying among 0.3, 0.6, and 0.9; the last one is the 30X simulated
sequencing of normal samples. The cells in the normal sample share the same
mutations, i.e., the germline mutations. 

For efficiency purposes, the tumour reads are not directly produced by rRACES.
Instead, the different coverages and purities are obtained by mixing two 200X
cohorts with purity 1 and 0. These two cohorts are called "tumour" and
"mixing_normal". Differently from the cohort normal sample, "mixing_normal"
simulates the reads of a sample whose cells contain all the tumour
pre-neoplastic mutations.

## Splitting "tumour" and "mixing_normal" into Lots

The `rRACES` functions `simulate_seq()` and `simulate_normal_seq()` can easily
produce the cohorts "tumour" and "mixing_normal". However, to speed up their
generation, we split each of them into 200 lots. Each of these lots simulates a
1X sequencing of the tumour at purity 1 (for "tumour") and purity 0 (for
"mixing_normal") using a different random seed so to guarantee the production
of a different set of reads.

The lot generation is distributed on a computational cluster managed by SLURM.
The cluster consists of 8 nodes equipped with an AMD EPYC 7H12 64-core
processor with 512GB of RAM and 352GB of local scratch.

## Producing a Set of Reads with a Specific Coverage and Purity

The number of reads $R$ to be simulated to reach a specific coverage $c$ is
$R = c\frac{G}{s}$ where $G$ and $s$ are the reference genome size and the read
size, respectively.

Some of the $R$ reads comes from "tumour" and some from "mixed_normal"
according to a probability that depends on the quantities of tumour and normal
DNA in the sample ($T$ and $N$, respectively). In particular, the probability
for one of the reads to come from a "mixed_normal" cell is $p_N=\frac{N}{T+N}$.

The `rRACES` method `PhylogeneticForest$get_samples_info()` returns the number
of tumour cells $t$ and the quantity of tumour DNA in the sample with respect
to a normal cell $D$, i.e., the number of normal cells $e$ needed to equal the
sample tumour DNA. Thus, 
$p_N = \frac{N}{T+N} = \frac{n D}{e D +n D} = \frac{n}{e+n}$ where $n$ is the
number of normal cells in the sample.

If $\pi$ is the aimed purity, then $n = t\frac{1-\pi}{\pi}$ is the number of
"mixing_normal" cells that should be added to the tumour sample to have purity
$\pi$. It follows that when the purity is $\pi$, the probability for one of the
reads to come from a "mixing_normal" cell is
$p_N=\frac{t(1-\pi)}{e\pi+t(1-\pi)}$.

Dealing with 200 uniformly distributed lots, we can sample the binomial
distribution $\mathcal{B}(p_N, R/200)$ 200 times and decide, lot by lot, how
many reads should be extracted from "mixing_normal" and, consequently, how many
from "tumour".

## Generating the Tumour Cohort

This directory contains 3 Python scripts that build the cohort of a SPN by
using its phylogenetic forest file on Orfeo.

To generate the cohort, copy all the Python scripts to a directory in your home
on Orfeo, log in on Orfeo, and execute the following command.

```{sh}
nohup python3 build_cohorts.py <SPN> <PHYLO_FOREST> <OUTPUT_DIR> -p EPYC -A cdslab &
```

where `<SPN>` is the name of the SPN, `<PHYLO_FOREST>` is the absolute path of
the phylogenetic forest, and `<OUTPUT_DIR>` is the absolute path of the output
directory.

For faster executions, both   `<PHYLO_FOREST>` and `<OUTPUT_DIR>` should be
located in one of the `/fast` subdirectories.

The Python scripts produce "tumour", "mixing_normal", and the normal sample
reads. Moreover, they also generate the FASTQ files and the corresponding
`sarek` files. 

> [!CAUTION]
> The content of the `/fast` subdirectories is deleted every 30 days
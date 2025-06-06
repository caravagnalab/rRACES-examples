#!/usr/bin/env Rscript
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

if (!require(sequenza)) stop("Package 'sequenza' missing\n.")

args <- commandArgs(TRUE)
input <- args[1]
output_prefix <- args[2]
gender <- args[3]

# as ascat
gam <- 70

# ploidy as ascat 
low_p <- 1.05
up_p <- 5.5

# complete spectrum of purity
high_ccf <- 1.0
low_ccf <- 0.1


if (gender == 'XX'){
  is_female = TRUE
} else {
  is_female = FALSE
}

output_dir = '.'
window = 1e5
overlap = 1
gamma = gam
kmin = 300
min_reads = 40
min_reads_normal = 10
min_reads_baf = 1
max_mut_types = 1
breaks = NULL
assembly = "hg19"
weighted_mean = TRUE
normalization_method = "mean"
segment_filter = 3e6
ratio_priority = FALSE
method = "baf"
low_cell = low_ccf
up_cell = high_ccf
low_ploidy = low_p
up_ploidy = up_p
CNt_max = 20

#  Define chromosomes to analyse (note these will subset to those that
# are available for sequenza:
chr_list <- c(1:22, "X")
if (gender != "XX") {
  chr_list <- c(chr_list, "Y")
}

chr_list <- c(chr_list, paste0("chr", chr_list))

# Extract sequenza data for model:
gc_stats <- gc.sample.stats(input)

modDat <- sequenza.extract(input,
                           window = window,
                           overlap = overlap,
                           gamma = gamma,
                           kmin = kmin,
                           min.reads = min_reads,
                           min.reads.normal = min_reads_normal,
                           min.reads.baf = min_reads_baf,
                           max.mut.types = max_mut_types,
                           chromosome.list = chr_list,
                           breaks = breaks,
                           assembly = assembly,
                           weighted.mean = weighted_mean,
                           normalization.method = normalization_method,
                           parallel = 8,
                           gc.stats = gc_stats
)

# Fit the model:
cells <- seq(low_ccf, high_ccf, 0.01)
plo <- seq(low_p, up_p, 0.1)
fit <- sequenza.fit(modDat,
                    female = is_female,
                    segment.filter = segment_filter,
                    cellularity = cells,
                    ploidy = plo,
                    XY = c(X = "chrX", Y = "chrY"),
                    ratio.priority = ratio_priority,
                    method = method
)

# Export the data:
sequenza.results(modDat,
                 fit,
                 output_prefix,
                 out.dir = output_dir,
                 female = is_female,
                 CNt.max = CNt_max,
                 XY = c(X = "chrX", Y = "chrY"),
                 ratio.priority = ratio_priority)
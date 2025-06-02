process SEQUENZAUTILS_RSEQZ {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-sequenza%3A3.0.0--r42h3342da4_5' :
        'biocontainers/r-sequenza%3A3.0.0--r42h3342da4_5' }"

    input:
    tuple val(meta), path(seqz_bin)

    output:
    tuple val(meta), path("*"), emit: rseqz
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def seqz_in = "${seqz_bin.toString().minus(".gz")}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${tissue}"
    
    """
    #!/usr/bin/env Rscript
    if (!require(sequenza)) stop("Package 'sequenza' missing\n.")
    
    input <- ${seqz_bin}
    output_prefix <- ${meta.id}
    gender <- ${meta.gender}
    ploidy <- params.ploidy
    ccf <- params.purity
    gam <- params.gamma
    
    if (ploidy == 7) {
      low_p <- 1
      up_p <- 7
    } else  {
      low_p <- ploidy - 0.5
      up_p <- ploidy + 0.5
    }
    
    if (ccf == 1) {
      high_ccf <- 1.0
      low_ccf <- 0.1
    } else {
      high_ccf <- ccf + 0.1
      low_ccf <- ccf - 0.1
    }
    
    if (gender == 'XX'){
      is_female = TRUE
    } else {
      is_female = FALSE
    }
    

    sequenzaAnalysis <- function(input,
                                 output_prefix,
                                 output_dir = '',
                                 window = 1e5,
                                 overlap = 1,
                                 gamma = gam,
                                 kmin = 300,
                                 min_reads = 40,
                                 min_reads_normal = 10,
                                 min_reads_baf = 1,
                                 max_mut_types = 1,
                                 breaks = NULL,
                                 assembly = "hg38",
                                 weighted_mean = TRUE,
                                 normalization_method = "mean",
                                 is_female = is_female,
                                 segment_filter = 3e6,
                                 ratio_priority = FALSE,
                                 method = "baf",
                                 low_cell = low_ccf,
                                 up_cell = high_ccf,
                                 low_ploidy = low_p,
                                 up_ploidy = up_p,
                                 CNt_max = 20) {
    
        #  Define chromosomes to analyse (note these will subset to those that
        # are available for sequenza:
        chr_list <- c(1:22, "X")
        if (gender != "XX") {
            chr_list <- c(chr_list, "Y")
        }
    
        chr_list <- c(chr_list, paste0("chr", chr_list))
    
        # Extract sequenza data for model:
        cat(sprintf("- Starting analysis for %s\n", input))
        cat("- Calculating gc-stats\n")
        gc_stats <- gc.sample.stats(input)
    
        cat("- Loading data\n")
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
        cat("- Fitting the model\n")
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
        cat("- Exporting results\n")
        sequenza.results(modDat,
            fit,
            output_prefix,
            out.dir = output_dir,
            female = is_female,
            CNt.max = CNt_max,
            XY = c(X = "chrX", Y = "chrY"),
            ratio.priority = ratio_priority
        )
    }
    
    sequenzaAnalysis(input, output_prefix)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def seqz_in = "${seqz_bin.toString().minus(".gz")}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${tissue}"
    """
    echo "analyse_cn_sequenza.R ${seqz_in} ${prefix} ${meta.id} ${gender} ${ploidy} ${purity} ${seq_gam}"
    mkdir ${tissue}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """
}

#' Convert Sequencing Results from Wide to Long Format
#'
#' This function takes sequencing results in wide format and converts them into a long format data frame.
#' It extracts sample names from column names, processes each sample separately, and then binds them together.
#' Finally, it renames and reorders columns to match the desired output format.
#'
#' @param seq_results A data frame containing sequencing results in wide format.
#' @return A data frame in long format with columns 'chr', 'from', 'to', 'ALT', 'NV', 'DP', 'VAF', and 'sample_name'.
#'
#' @examples
#' # Example data frame in wide format
#' seq_results <- data.frame(chr = c("chr1", "chr2"),
#'                            chr_pos = c(100, 200),
#'                            ref = c("A", "C"),
#'                            alt = c("T", "G"),
#'                            sample1.VAF = c(0.1, 0.2),
#'                            sample2.VAF = c(0.3, 0.4))
#' # Convert to long format
#' seq_to_long(seq_results)
#'
#' @export
seq_to_long <- function(seq_results) {
  # Extract sample names from column names
  sample_names <- strsplit(colnames(seq_results)[grepl(".VAF", colnames(seq_results), fixed = TRUE)], ".VAF") %>% unlist()

  # Process each sample separately to create a list of data frames
  seq_df <- lapply(sample_names, function(sn) {
    # Select relevant columns for the current sample
    cc <- c("chr", "chr_pos", "ref", "alt", "causes", "classes", colnames(seq_results)[grepl(paste0(sn, "."), colnames(seq_results), fixed = TRUE)])

    # Rename columns and add sample_name column
    seq_results[, cc] %>%
      `colnames<-`(c("chr", "chr_pos", "ref", "alt", "causes", "classes", "occurences", "coverage", "VAF")) %>%
      dplyr::mutate(sample_name = sn)
  }) %>% do.call("bind_rows", .)

  # Rename and reorder columns
  seq_df %>%
    dplyr::rename(chr = chr, from = chr_pos, DP = coverage, NV = occurences, ALT = alt) %>%
    dplyr::mutate(to = from)
}

#' Plot Depth Ratio (DR) Genome-wide
#'
#' This function generates a plot showing the DEpth Ratio (DR) genome-wide for a specific sample.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot. Default is all autosomes and sex chromosomes.
#' @param n The number of mutations to sample for plotting. Default is 5000.
#' @return A ggplot2 object showing the DR distribution across the genome.
#'
#' @export
plot_DR_gw <- function(seq_res, sample, chromosomes = paste0(c(1:22, "X", "Y")), n = 5000) {
  if (any(!chromosomes %in% paste0(c(1:22, "X", "Y"))))
    stop("Invalid chromosome passed as input")

  data <- seq_to_long(seq_res) %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(classes == "germinal") %>%
    dplyr::mutate(mut_id = paste(chr, from, to, sep = ":"))

  normal_data <- data %>% dplyr::filter(sample_name == "normal_sample")
  if (!sample %in% unique(data$sample_name)) {
    stop(paste("Inserted 'sample' is not available, available samples are:",
               paste0(unique(data$sample_name), collapse = ", ")))
  }
  tumour_data <- data %>% dplyr::filter(sample_name == sample)

  d <- tumour_data %>%
    dplyr::left_join(normal_data, suffix = c(".tumour", ".normal"), by = "mut_id") %>%
    dplyr::mutate(DR = DP.tumour / DP.normal) %>%
    dplyr::group_by(chr.tumour) %>%
    dplyr::sample_n(n) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(chr.tumour), desc(from.tumour)) %>%
    dplyr::mutate(abs_pos = 1:n())

  chr_limits <- d %>%
    dplyr::group_by(chr.tumour) %>%
    dplyr::filter(abs_pos == min(abs_pos)) %>%
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(d$abs_pos))

  chr_means <- d %>%
    dplyr::group_by(chr.tumour) %>%
    dplyr::summarise(mean_pos = (max(abs_pos) + min(abs_pos)) / 2) %>%
    dplyr::pull(mean_pos)

  d %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = abs_pos, y = DR)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid") +
    ggplot2::geom_point(size = 0.2, alpha = 0.2) +
    ggplot2::labs(x = "", y = "DR") +
    ggplot2::lims(y = c(0, NA)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr.tumour)) +
    ggplot2::geom_hline(yintercept = stats::median(d$DR), col = "indianred", linetype = "dashed")
}


#' Plot B-Allele Frequency (BAF) Genome-wide
#'
#' This function generates a plot showing the B-Allele Frequency (BAF) genome-wide for a specific sample.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot.
#'                    Default is all autosomes and sex chromosomes.
#' @param cuts A numeric vector specifying the range of BAF values to include in the plot.
#'             Default is c(0.01, 0.99).
#' @param n The number of mutations to sample for plotting. Default is 5000.
#' @return A ggplot2 object showing the BAF distribution across the genome.
#'
#' @export
plot_BAF_gw <- function(seq_res, sample, chromosomes = paste0(c(1:22, "X", "Y")), cuts = c(0.01, 0.99), n = 5000) {
  if (any(!chromosomes %in% paste0(c(1:22, "X", "Y"))))
    stop("Invalid chromosome passed as input")

  data <- seq_to_long(seq_res) %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(classes == "germinal")

  if (!sample %in% unique(data$sample_name)) {
    stop(paste("Inserted 'sample' is not available, available samples are:",
               paste0(unique(data$sample_name), collapse = ", ")))
  }
  tumour_data <- data %>% dplyr::filter(sample_name == sample)

  d <- tumour_data %>%
    dplyr::group_by(chr) %>%
    dplyr::sample_n(n) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(chr), desc(from)) %>%
    dplyr::mutate(abs_pos = 1:n()) %>%
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))

  chr_limits <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(abs_pos == min(abs_pos)) %>%
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(d$abs_pos))

  chr_means <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(mean_pos = (max(abs_pos) + min(abs_pos)) / 2) %>%
    dplyr::pull(mean_pos)

  d %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = abs_pos, y = VAF)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid") +
    ggplot2::geom_point(size = 0.2, alpha = 0.2) +
    ggplot2::labs(x = "", y = "BAF") +
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr)) +
    ggplot2::geom_hline(yintercept = stats::median(d$VAF), col = "indianred", linetype = "dashed")
}


#' Plot Variant Allele Frequency (VAF) Genome-wide
#'
#' This function generates a plot showing the Variant Allele Frequency (VAF) genome-wide for a specific sample.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot.
#'                    Default is all autosomes and sex chromosomes.
#' @param cuts A numeric vector specifying the range of VAF values to include in the plot.
#'             Default is c(0.05, 0.99).
#' @param n The number of mutations to sample for plotting. Default is 1000.
#' @return A ggplot2 object showing the VAF distribution across the genome.
#'
#' @export
plot_VAF_gw <- function(seq_res, sample, chromosomes = paste0(c(1:22, "X", "Y")), cuts = c(0.05, 0.99), n = 1000) {
  if (any(!chromosomes %in% paste0(c(1:22, "X", "Y"))))
    stop("Invalid chromosome passed as input")

  data <- seq_to_long(seq_res) %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(classes != "germinal")

  if (!sample %in% unique(data$sample_name)) {
    stop(paste("Inserted 'sample' is not available, available samples are:",
               paste0(unique(data$sample_name), collapse = ", ")))
  }
  tumour_data <- data %>% dplyr::filter(sample_name == sample)

  d <- tumour_data %>%
    dplyr::arrange(desc(chr), desc(from)) %>%
    dplyr::group_by(chr) %>%
    dplyr::sample_n(n, replace = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(abs_pos = 1:n()) %>%
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))

  chr_limits <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(abs_pos == min(abs_pos)) %>%
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(d$abs_pos))

  chr_means <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(mean_pos = (max(abs_pos) + min(abs_pos)) / 2) %>%
    dplyr::pull(mean_pos)

  d %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = abs_pos, y = VAF)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid") +
    ggplot2::geom_point(size = 0.5, alpha = 0.4) +
    ggplot2::labs(x = "", y = "VAF") +
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr)) +
    ggplot2::geom_hline(yintercept = stats::median(d$VAF), col = "indianred", linetype = "dashed")
}

#' Plot Histogram of Variant Allele Frequency (VAF)
#'
#' This function generates a histogram showing the distribution of Variant Allele Frequency (VAF) across samples and chromosomes.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot.
#'                    Default is all autosomes and sex chromosomes.
#' @param samples A character vector specifying the sample names to include in the plot.
#'                Default is NULL, which includes all samples except the "normal_sample".
#' @param colour_by A character indicating whether to color the histogram bars by "causes" or "classes".
#'                  Default is "causes".
#' @param cuts A numeric vector specifying the range of VAF values to include in the plot.
#'             Default is c(0, 1).
#' @return A ggplot2 object showing the histogram of VAF.
#'
#' @export
plot_histogram_vaf <- function(
    seq_res,
    chromosomes = paste0(c(1:22, "X", "Y")),
    samples = NULL,
    colour_by = "causes",
    cuts = c(0, 1)
) {
  if (!(colour_by %in% c("classes", "causes")))
    stop("Colour_by parameter can be either 'causes' or 'classes'")
  if (any(!chromosomes %in% paste0(c(1:22, "X", "Y"))))
    stop("Invalid chromosomes passed as input")

  data <- seq_to_long(seq_res) %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(classes != "germinal") %>%
    dplyr::filter(sample_name != "normal_sample") %>%
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))

  if (!is.null(samples)) {
    if (any(!samples %in% unique(data$sample_name)))
      stop("Invalid sample name in samples parameter")
    data <- data %>% dplyr::filter(sample_name %in% samples)
  }

  data$col <- data[[colour_by]]

  data %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = VAF, col = col, fill = col)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::xlim(x = c(-0.01, 1.01)) +
    ggplot2::facet_grid(sample_name ~ chr, scales = 'free') +
    ggplot2::theme_bw() +
    ggplot2::labs(col = colour_by, fill = colour_by) +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot Marginals of Variant Allele Frequency (VAF)
#'
#' This function generates scatter plots showing the marginal distributions of Variant Allele Frequency (VAF) for pairs of samples on a specific chromosome.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param chromosome A character specifying the chromosome to analyze.
#' @param colour_by A character indicating whether to color the scatter points by "causes" or "classes".
#'                  Default is "causes".
#' @param samples A character vector specifying the sample names to include in the plot.
#'                Default is NULL, which includes all samples except the "normal_sample".
#' @param cuts A numeric vector specifying the range of VAF values to include in the plot.
#'             Default is c(0, 1).
#' @return A list of ggplot2 objects showing scatter plots of VAF marginals for pairs of samples.
#'
#' @export
plot_marginals <- function(seq_res, chromosome, colour_by = "causes", samples = NULL, cuts = c(0, 1)) {
  if (!(colour_by %in% c("classes", "causes")))
    stop("Colour_by parameter can be either 'causes' or 'classes'")
  if (!chromosome %in% paste0(c(1:22, "X", "Y")))
    stop("Invalid chr passed as input")

  data <- seq_to_long(seq_res) %>%
    dplyr::filter(chr == chromosome) %>%
    dplyr::filter(classes != "germinal") %>%
    dplyr::filter(sample_name != "normal_sample") %>%
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))

  if (!is.null(samples)) {
    if (any(!samples %in% unique(data$sample_name)))
      stop("Invalid sample name in samples parameter")
    data <- data %>% dplyr::filter(sample_name %in% samples)
  }

  combinations <- utils::combn(unique(data$sample_name), m = 2)

  lapply(1:ncol(combinations), function(i) {
    couple <- combinations[, i]
    d1 <- data %>% dplyr::filter(sample_name == couple[1]) %>% dplyr::mutate(mut_id = paste(chr, from, to, sep = ":"))
    d2 <- data %>% dplyr::filter(sample_name == couple[2]) %>% dplyr::mutate(mut_id = paste(chr, from, to, sep = ":"))

    djoin <- dplyr::full_join(d1, d2, by = 'mut_id')
    djoin$col <- djoin[[paste0(colour_by, ".x")]]

    djoin %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = VAF.x, y = VAF.y, col = col)) +
      ggplot2::geom_point() +
      ggplot2::xlim(c(-0.01, 1.01)) +
      ggplot2::ylim(c(-0.01, 1.01)) +
      ggplot2::labs(x = couple[1], y = couple[2], col = colour_by) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")
  })
}


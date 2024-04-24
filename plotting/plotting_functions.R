
seq_to_long <- function(seq_results) {
  # get names of samples
  sample_names <- strsplit(colnames(seq_results)[grepl(".VAF", colnames(seq_results), fixed = TRUE)], ".VAF") %>% unlist()
  
  seq_df <- lapply(sample_names, function(sn) {
    cc <- c("chr", "chr_pos", "ref", "alt", "causes", "classes", colnames(seq_results)[grepl(paste0(sn, "."), colnames(seq_results), fixed = T)])
    
    seq_results[, cc] %>%
      `colnames<-`(c("chr", "chr_pos", "ref", "alt", "causes", "classes", "occurences", "coverage", "VAF")) %>%
      dplyr::mutate(sample_name = sn)
  }) %>% do.call("bind_rows", .)
  
  seq_df %>%
    dplyr::rename(chr = chr, from = chr_pos, DP = coverage, NV = occurences, ALT = alt) %>%
    dplyr::mutate(to = from)
  #dplyr::select(chr, from, to, ALT, NV, DP, VAF, sample_name)
}

# plot gw stuff ####
plot_DR_gw <- function(seq_res, sample, chromosomes = paste0(c(1:22, "X", "Y"))) {
  data <- seq_to_long(seq_res) %>% 
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>% 
    dplyr::filter(chr %in% chromosomes) %>% 
    dplyr::filter(classes == "germinal")
    
  normal_data <- data %>% dplyr::filter(sample_name == "normal_sample")
  tumour_data <- data %>% dplyr::filter(sample_name == sample)  
  
  d <- tumour_data %>%
    dplyr::left_join(normal_data, suffix = c(".tumour", ".normal"), by="from") %>%
    dplyr::mutate(DR = DP.tumour / DP.normal) %>% 
    dplyr::arrange(desc(chr.tumour), desc(from)) %>% 
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
    ggplot2::ggplot(mapping = ggplot2::aes(x=abs_pos, y=DR)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid") +
    ggplot2::geom_point(size = .2, alpha = .2) +
    ggplot2::labs(x = "", y = "DR") +
    ggplot2::lims(y = c(0,NA)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr.tumour)) +
    ggplot2::geom_hline(yintercept = stats::median(d$DR), col = "indianred", linetype = "dashed")
}

plot_BAF_gw <- function(seq_res, sample, chromosomes = paste0(c(1:22, "X", "Y")), cuts = c(.01, .99)) {
  data <- seq_to_long(seq_res) %>% 
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>% 
    dplyr::filter(chr %in% chromosomes) %>% 
    dplyr::filter(classes == "germinal")
  
  tumour_data <- data %>% dplyr::filter(sample_name == sample)  
  
  d <- tumour_data %>%
    dplyr::arrange(desc(chr), desc(from)) %>%
    dplyr::mutate(abs_pos = 1:n()) %>% 
    dplyr::filter()
  
  chr_limits <- d %>% 
    dplyr::group_by(chr) %>% 
    dplyr::filter(abs_pos == min(abs_pos)) %>% 
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(d$abs_pos))
  
  chr_means <- d %>% 
    dplyr::group_by(chr) %>% 
    dplyr::summarise(mean_pos = (max(abs_pos) + min(abs_pos)) / 2) %>% 
    dplyr::pull(mean_pos)
  
  d <- d %>% dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))
  
  d %>% 
    ggplot2::ggplot(mapping = ggplot2::aes(x=abs_pos, y=VAF)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid") +
    ggplot2::geom_point(size = .2, alpha = .2) +
    ggplot2::labs(x = "", y = "BAF") +
    ggplot2::lims(y = c(0,1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr.tumour)) +
    ggplot2::geom_hline(yintercept = stats::median(d$VAF), col = "indianred", linetype = "dashed")
}

plot_VAF_gw <- function(seq_res, sample, chromosomes = paste0(c(1:22, "X", "Y")), cuts = c(.05, .99)) {
  data <- seq_to_long(seq_res) %>% 
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>% 
    dplyr::filter(chr %in% chromosomes) %>% 
    dplyr::filter(classes != "germinal")
  
  tumour_data <- data %>% dplyr::filter(sample_name == sample)  
  
  d <- tumour_data %>%
    dplyr::arrange(desc(chr), desc(from)) %>%
    dplyr::mutate(abs_pos = 1:n()) %>% 
    dplyr::filter()
  
  chr_limits <- d %>% 
    dplyr::group_by(chr) %>% 
    dplyr::filter(abs_pos == min(abs_pos)) %>% 
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(d$abs_pos))
  
  chr_means <- d %>% 
    dplyr::group_by(chr) %>% 
    dplyr::summarise(mean_pos = (max(abs_pos) + min(abs_pos)) / 2) %>% 
    dplyr::pull(mean_pos)
  
  d <- d %>% dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))
  
  d %>% 
    ggplot2::ggplot(mapping = ggplot2::aes(x=abs_pos, y=VAF)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid") +
    ggplot2::geom_point(size = .5, alpha = .4) +
    ggplot2::labs(x = "", y = "VAF") +
    ggplot2::lims(y = c(0,1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr.tumour)) +
    ggplot2::geom_hline(yintercept = stats::median(d$VAF), col = "indianred", linetype = "dashed")
}

# Plot data histogram ####
plot_histogram_vaf <- function(
    seq_res, 
    chromosomes = paste0(c(1:22, "X", "Y")), 
    samples = NULL, 
    colour_by = "causes",
    cuts = c(0,1)
) {
  if (!(colour_by %in% c("classes", "causes"))) stop("colour_by parameter can be either 'causes' or 'classes'")
  if (any(!chromosomes %in% paste0(c(1:22, "X", "Y")))) stop("invalid chromosomes passed as input")
  
  data <- seq_to_long(seq_res) %>% 
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>% 
    dplyr::filter(chr %in% chromosomes) %>% 
    dplyr::filter(classes != "germinal") %>% 
    dplyr::filter(sample_name != "normal_sample") %>% 
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))
  
  if (!is.null(samples)) { 
    if (any(!samples %in% unique(data$sample_name))) stop("invalid sample name in samples parameter")  
    data <- data %>% dplyr::filter(sample_name %in% samples) 
  }
  
  data$col = data[[colour_by]]
  
  data %>% 
    ggplot2::ggplot(mapping = ggplot2::aes(x=VAF, col=col, fill=col)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::xlim(x=c(-0.01, 1.01)) +
    ggplot2::facet_grid(sample_name ~ chr, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::labs(col=colour_by, fill=colour_by) +
    ggplot2::theme(legend.position = "bottom")  
}

# plot marginals ####
plot_marginals <- function(seq_res, chromosome, colour_by="causes", samples=NULL) {
  if (!(colour_by %in% c("classes", "causes"))) stop("colour_by parameter can be either 'causes' or 'classes'")
  if (!chromosome %in% paste0(c(1:22, "X", "Y"))) stop("invalid chr passed as input")
  
  data <- seq_to_long(seq_res) %>% 
    dplyr::filter(chr == chromosome) %>% 
    dplyr::filter(classes != "germinal") %>% 
    dplyr::filter(sample_name != "normal_sample") %>% 
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))
  
  if (!is.null(samples)) { 
    if (any(!samples %in% unique(data$sample_name))) stop("invalid sample name in samples parameter")  
    data <- data %>% dplyr::filter(sample_name %in% samples) 
  }
  
  combinations <- utils::combn(unique(data$sample_name), m = 2)
  
  lapply(1:ncol(combinations), function(i) {
    couple <- combinations[,i]
    d1 <- data %>% dplyr::filter(sample_name == couple[1]) %>% dplyr::mutate(mut_id = paste(chr, from, to, sep = ":"))
    d2 <- data %>% dplyr::filter(sample_name == couple[2]) %>% dplyr::mutate(mut_id = paste(chr, from, to, sep = ":"))
    
    djoin <- dplyr::full_join(d1, d2, by = 'mut_id')
    djoin$col = djoin[[paste0(colour_by, ".x")]]
    
    djoin %>% 
      ggplot2::ggplot(mapping = ggplot2::aes(x = VAF.x, y = VAF.y, col = col)) +
      ggplot2::geom_point() +
      ggplot2::xlim(c(-0.01,1.01)) +
      ggplot2::ylim(c(-0.01,1.01)) +
      ggplot2::labs(x = couple[1], y=couple[2], col=colour_by) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")  
  })
}

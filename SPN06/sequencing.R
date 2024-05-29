
# load phylogenetic forest
phylo_forest <- load_phylogenetic_forest("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/phylo_forest.sff")

seq_results <- simulate_seq(phylo_forest, coverage = 2.5)


seq_results <- simulate_seq(phylo_forest, coverage = 2.5, write_SAM = FALSE)


data <- seq_to_long(seq_results)

BASE_FACTOR <- 1e6



data %>% dplyr::filter(chr == "22") %>% pull(sample_name) %>% table()
curr_chr <- "22"
for (curr_chr in data$chr %>% unique()) {
  normal_data <- data %>%
    dplyr::filter(chr == curr_chr) %>%
    dplyr::filter(sample_name == "normal_sample")
  
  d <- data %>%
    dplyr::filter(chr == curr_chr) %>%
    dplyr::filter(sample_name != "normal_sample")
  
  dd <- d %>%
    dplyr::left_join(normal_data, suffix = c(".tumour", ".normal"), by="from") %>%
    dplyr::mutate(DR = DP.tumour / DP.normal)
  
  DR_plot <- dd %>%
    dplyr::filter(classes.normal == "germinal") %>%
    ggplot(mapping = aes(x = from / BASE_FACTOR, y=DR)) +
    geom_point(size = .2, alpha = .2) +
    labs(x = "", y = "DR") +
    ggtitle(paste0("Chromsome ", curr_chr)) +
    facet_wrap(~ sample_name.tumour) +
    theme_bw()
  
  BAF_plot <- dd %>%
    dplyr::filter(classes.normal == "germinal", VAF.tumour < 0.98) %>%
    ggplot(mapping = aes(x = from / BASE_FACTOR, y=VAF.tumour)) +
    geom_point(size = .2, alpha = .2) +
    labs(x = "", y = "BAF") +
    lims(y = c(0,1)) +
    facet_wrap(~ sample_name.tumour) +
    theme_bw()
  
  VAF_plot <- dd %>%
    dplyr::filter(classes.normal != "germinal") %>%
    ggplot(mapping = aes(x = from / BASE_FACTOR, y=VAF.tumour)) +
    geom_point(size = .2, alpha = .2) +
    labs(x = "MegaBase", y = "VAF") +
    lims(y = c(0,1)) +
    facet_wrap(~ sample_name.tumour) +
    theme_bw()
  
  VAF_spectrum <- dd %>%
    dplyr::filter(classes.normal != "germinal", VAF.tumour > 0.02) %>%
    ggplot(mapping = aes(x = VAF.tumour)) +
    geom_histogram(bins = 100) +
    labs(x = "VAF", y = "Counts") +
    lims(x = c(0,1)) +
    facet_wrap(~ sample_name.tumour) +
    theme_bw()
  
  pp <- DR_plot / BAF_plot / VAF_plot / VAF_spectrum
  ggsave(paste0("/Users/azadsadr/Documents/packages/rRACES-examples/SPN06/chromosomes/", curr_chr, ".png"), dpi = 400, width = 6, height = 8, plot = pp)
}





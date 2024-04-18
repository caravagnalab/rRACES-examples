rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils.R")

data <- readRDS("data/sequencing.rds")
BASE_FACTOR <- 1e6

data <- seq_to_long(data)
nrow(data) / 1e6

data %>% dplyr::filter(chr == "1") %>% pull(sample_name) %>% table()
curr_chr <- "1"
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
    geom_point(size = .5) +
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
  ggsave(paste0("tissue/chromosomes/", curr_chr, ".png"), dpi = 400, width = 6, height = 8, plot = pp)
}


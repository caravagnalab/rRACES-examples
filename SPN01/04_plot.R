rm(list = ls())
library(rRACES)
library(dplyr)
library(ggplot2)
library(patchwork)


### ALLELE-SPEC ###
phylo_forest <- load_phylogenetic_forest("data/phylo_forest_WGD.sff")
samples_name <- phylo_forest$get_samples_info()[["name"]]

plot_ASCNP <- function(bulk_allelic_fragm, sample){
  d <- as_tibble(bulk_allelic_fragm) %>%
    dplyr::mutate(chr = as.numeric(chr)) %>%
    dplyr::arrange(chr, begin, end)

  max_pos <- aggregate(end ~ chr, data = d, max)
  # Step 2: Calculate the cumulative sum of chromosome lengths
  max_pos$cumulative_length <- c(0, cumsum(as.numeric(max_pos$end)[-length(max_pos$end)]))

  # Step 3: Merge cumulative lengths back with the original data
  d <- merge(d, max_pos[, c("chr", "cumulative_length")], by = "chr")

  # Step 4: Convert to absolute coordinates
  d$abs_begin <- d$begin + d$cumulative_length
  d$abs_end <- d$end + d$cumulative_length

  # View the result
  d <- d[order(d$chr, d$begin), ]  # Order by chromosome and begin position
  d_long <- tidyr::gather(d, key = "allele_type", value = "count", major, minor)

  chr_limits <- d_long %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(abs_begin == min(abs_begin)) %>%
    dplyr::pull(abs_begin)
  chr_limits <- c(chr_limits, max(d_long$abs_begin))

  chr_means <- d_long %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(mean_pos = (max(abs_end) + min(abs_begin)) / 2) %>%
    dplyr::pull(mean_pos)

  aspc_plot <- d_long %>%
    ggplot(aes(x = abs_begin, xend = abs_end, y = count, yend = count, color = allele_type)) +
    geom_segment(size = 1, position = position_dodge(width = 0.2)) +
    scale_color_manual(values = c("major" = "red", "minor" = "blue")) +
    labs(
      title = sample,
      y = "Allele Count",
      color = "Allele Type"
    ) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid",
                        alpha = 0.4) +
    ggplot2::geom_point(size = 0.2, alpha = 0.3) +
    ggplot2::labs(x = "", y = "Allele") +
    ggplot2::lims(y = c(0, NA)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means,
                                labels = unique(d$chr)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 50, by = 1))
  return(aspc_plot)
}

#seq_res <- readRDS("data/seq_results_80X_with_error_paired.rds")
seq_res <- readRDS("data/seq_results_WGD_muCNA_5e11.rds")
### VAFS ###
#samples <- c('SPN01_Sample_1','SPN01_Sample_2','SPN01_Sample_3')
samples <- c("Sample_A","Sample_B","Sample_C")
g_seq <- seq_res %>% filter(classes=="germinal")
for (sample in samples) {
  print(sample)
  #pdf(paste0(sample, '_report.pdf'), width = 8, height = 10)
  #pdf(paste0("tissue/chromosomes/", sample, '_report.pdf'), width = 8, height = 10)
  p1 <- plot_VAF(g_seq, sample)
  p2 <- plot_BAF(g_seq, sample)
  p3 <- plot_DR(g_seq, sample)
  bulk_sample <- phylo_forest$get_bulk_allelic_fragmentation(sample)
  p4<- plot_ASCNP(bulk_allelic_fragm = bulk_sample, sample = sample)

  p <- p1 / p2 / p3 / p4
  ggsave(paste0("plots/",sample, '_report_WGD.png'),plot = p, dpi=300, width = 8, height = 14)
  #p %>% print()
  #dev.off()
}

### MARGINALS ###
pdf("plots/chromosome_vaf_marginals_report_WGD.pdf", width = 16, height = 5)
seq_res$chr %>% unique()
s_seq <- seq_res %>% filter(classes!="germinal")
for (c in unique(seq_res$chr)) {
    print(c)
    p_marg <- plot_VAF_marginals(s_seq, chromosomes = c, samples = samples, labels = s_seq["classes"])
    p_hist <- plot_VAF_histogram(s_seq, chromosomes = c, samples = samples, labels = s_seq["classes"], cuts = c(0.02, 1))
    p <- (wrap_plots(p_marg, ncol=4)+p_hist)+ plot_layout(guides = 'collect') & theme(legend.position = 'bottom') & ggtitle(paste("Chromosome", c))
    #p_marg <- lapply(p_marg, function(p) p + ggtitle(paste("Chromosome", c)))
    #p <- wrap_plots(list(p_marg,p_hist),ncol = 3, nrow=2) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
    print(p)
}

dev.off()
#marginals <- lapply(unique(seq_res$chr), function(c) {
#			    print(c)
#			    p_marg <-plot_VAF_marginals(s_seq, chromosomes=c, samples = samples, labels = s_seq["classes"])
#			    p_marg <- lapply(p_marg, function(p) p + ggtitle(paste("Chromosome", c)))
#			    p <- wrap_plots(p_marg, ncol = 3) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
#			    print(p)
#			    #ggsave(paste0("chr_",c,"_vaf_marginals.pdf"), dpi=300, width = 16, height = 8,plot=p)
#})


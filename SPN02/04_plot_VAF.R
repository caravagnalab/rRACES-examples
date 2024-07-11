library(rRACES)
library(tidyverse)
library(patchwork)

seq_res = readRDS("SPN02_seq_80x.rds")
# seq_res = seq_to_long(seq_res)

# removed germline mutations
# seq_res_somatic <- seq_res %>% 
#     filter(classes != "germinal")

seq_res = bind_rows(seq_res)
samples = c("A", "B")
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

seq_res_filt = seq_res %>% filter(classes != "germinal")
seq_res_germ = seq_res %>% filter(classes == "germinal")
lapply(samples, function(s) {
    vaf = plot_VAF(seq_res_filt, s)
    baf = plot_BAF(seq_res_germ, s)
    dr = plot_DR(seq_res_germ, s)

    p = vaf / baf / dr

    ggsave(paste0("plot/", s, "_report.png"), plot = p, width = 8, height = 10)
}) 

pdf("plot/chromosome_vaf_marginals_report.pdf", width = 16, height = 5)
# seq_res_filt$chr %>% unique()
# s_seq <- seq_res %>% filter(classes!="germinal")
for (c in unique(seq_res_filt$chr)) {
    print(c)
    p_marg <- plot_VAF_marginals(seq_res_filt, chromosomes = c, samples = samples, labels = seq_res_filt["classes"])
    p_hist <- plot_VAF_histogram(seq_res_filt, chromosomes = c, samples = samples, labels = seq_res_filt["classes"], cuts = c(0.02, 1))
    p <- (wrap_plots(p_marg, ncol=2)+p_hist)+ plot_layout(guides = 'collect') & theme(legend.position = 'bottom') & ggtitle(paste("Chromosome", c))
    #p_marg <- lapply(p_marg, function(p) p + ggtitle(paste("Chromosome", c)))
    #p <- wrap_plots(list(p_marg,p_hist),ncol = 3, nrow=2) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
    print(p)
}

dev.off()
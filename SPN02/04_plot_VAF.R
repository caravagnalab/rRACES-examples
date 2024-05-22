library(rRACES)
library(tidyverse)
library(patchwork)

seq_res = readRDS("test/SPN02/sequencing_wg_80x.rds")
# seq_res = seq_to_long(seq_res)

# removed germline mutations
# seq_res_somatic <- seq_res %>% 
#     filter(classes != "germinal")

seq_res = bind_rows(seq_res)
samples = c("A", "B")
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

# plot VAF one chromosome at time
vaf_gw = lapply(samples, function(x) {
    vaf_gw = lapply(chromosomes, function(c) {
        plot_VAF_gw(seq_res, 
            sample = x, 
            chromosomes = c(c))
    })
    # vaf_gw_plot = wrap_plots(vaf_gw, ncol = 1, nrow = length(vaf_gw))
    # pdf(paste0("test/SPN02/vaf_gw_", x, ".pdf"))
    # vaf_gw
    # dev.off()
    # ggsave(filename = paste0("test/SPN02/vaf_gw_", x, ".pdf"), plot = vaf_gw_plot)
})
names(vaf_gw) = samples

pdf("test/SPN02/vaf_gw_a.pdf")
vaf_gw[[1]]
dev.off()
pdf("test/SPN02/vaf_gw_b.pdf")
vaf_gw[[2]]
dev.off()

# lapply(samples, function(x) {
plot = plot_histogram_vaf(seq_res, 
    # sample = x, 
    colour_by = "causes", 
    cuts = c(-0.1, 1.01)
    )
ggsave(filename = paste0("test/SPN02/vaf_hist.pdf"), plot = plot, width = 20)
# }) 
# make it bigger

marginals = lapply(chromosomes, function(x) {
    plot_marginals(
        seq_res, 
        chromosome = c(x), 
        cuts = c(-0.1, 1.01))
})

pdf("test/SPN02/marginal_vaf.pdf")
marginals
dev.off()
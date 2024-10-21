rm(list = ls())
library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CNAqc)
#source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN01/scripts/my_functions/plot_genome_wide.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/plotting/spn_blueprint/utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/plotting/spn_blueprint/plot_genome_wide.R")
seq_results <- readRDS("new_simulation/sim5/monday_meeting/WGD_post/seq_res_new.rds")
phylo_forest <- load_phylogenetic_forest("new_simulation/sim5/monday_meeting/WGD_post/phylo_forest.sff")

samples <- phylo_forest$get_samples_info()[["name"]] %>% sort()
print(samples)
pdf("genome_wide.pdf", width = 10, height = 12)
plots_gw <- list()

cov = 80
error_rate = 1e-3
tumour_type = "COAD"
germline_sub = "default"
#rownames <- c("row1", "row2", "row3")
#colnames <- c("col1", "col2", "col3")
#
## Create an empty dataframe
#df <- data.frame(matrix(ncol = length(colnames), nrow = length(rownames)))

#plots_insp <- list()
for (i in samples){
	x <- races2cnaqc(seq_results=seq_results,
			 phylo_forest=phylo_forest,
			 sample_id=i,
			 ref="GRCh38",
			 purity=1)
	print(x$cna)
	saveRDS(x,file=paste0("CNAqc_",i,".rds"))
        #g_seq <- seq_results %>% filter(classes=="germinal")
        #seg <- plot_segments(x = x,highlight = c("1:0","1:1","2:1","2:2"))
        #dr <- plot_DR_n(seq_res = g_seq,sample = i)
        #baf <- plot_BAF_n(seq_res = g_seq,sample = i)
        #vaf_histo <- plot_data_histogram(x, which = 'VAF')
        #dp_histo <- plot_data_histogram(x, which = 'DP')

        #p1 <- vaf_histo/dp_histo+plot_layout(design="AAAAA\nAAAAA\nAAAAA\nBBBBB\nBBBBB\nBBBBB")
        #p2 <- dr/baf/seg+plot_layout(design="AAAAA\nAAAAA\nBBBBB\nBBBBB\nCCCCC\nCCCCC")


	gw <- genome_wide_plots(x,seq_results,i)
	plots_gw[[i]] <- wrap_plots(gw[[1]],gw[[2]])+
		plot_layout(design="AAAAA\nAAAAA\nBBBBB\nBBBBB")+
		plot_annotation(title = i, subtitle = paste0("Simulated coverage: ", cov,
							     "\nSimulated purity: ", x$purity,
							     "\nSequencing Errore rate: ", error_rate,
							     "\nTumour type: ", tumour_type,
							     "\nGermline subject :", germline_sub))
	#plots_insp[[i]] <-inspect_segment(x)+ggtitle(label=i)
	print(plots_gw[[i]])
	seg <- wrap_plots(gw[[3]])
	print(seg)
}





#saveRDS(plots_gw,"page3_gw_1.rds")
#plot_vafs <- wrap_plots(plots_insp,ncol =1,nrow=length(samples))

#saveRDS(plot_vafs,"page3_vafs.rds")
dev.off()


chroms <- c("22","5","17","3")
#chroms <- unique(seq_results$chr)
#chroms <- c("1","2","3","4","5","6","7","8", "9","10","11","12","13","14","15",
#	    "16","17","18","19","20","21","22","X","Y")
pdf("chromosome_vafs_1.pdf",width = 8, height = 10)
vaf_plots <- list()
for (c in chroms){
        vaf_plots[[c]] <- squareplot(seq_res=seq_results,
			       samples_list=samples,
			       chrom=c)
	print(c)
	print(vaf_plots[[c]])
	#print(c)
}
dev.off()
#saveRDS(vaf_plots,"page4_vafs.rds")

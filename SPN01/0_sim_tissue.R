rm(list = ls())
library(rRACES)
library(dplyr)
library(patchwork)
library(ggplot2)
library(english)
library(cli)
library(ggrepel)
library(ggpubr)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/plotting/spn_blueprint/plot_tissue_dynamics.R")
#source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/rRACES-examples/SPN02/plotting_sampling.R")

dir <- getwd()
set.seed(12345)
sim <- SpatialSimulation(name = 'SPN01', seed = 12345)
tissue_size <- sim$get_tissue_size()
sim$history_delta <- 1
sim$death_activation_level <- 50

tissue_plots <- list()
state_plots <- list()

clone_rates <- c(0.08,0.3,1,3)
names(clone_rates) <- c("Clone 1","Clone 2","Clone 3","Clone 4")

time_intervals <- 50


# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = 0.08, death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 500)

t1 <- plot_tissue(sim)+annotate(geom="text",x=900,y=900,label=paste0("Time: ",as.character(round(sim$get_clock(),2))))
s1 <- plot_state(sim)

# Clone 2
sim$add_mutant(name = "Clone 2", growth_rates = 0.3, death_rates = 0.01)
sim$update_rates(species = "Clone 1", rates = c(growth = 0.06, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size(species = 'Clone 2', 700)

t2 <- plot_tissue(sim)+annotate(geom="text",x=900,y=900,label=paste0("Time: ",as.character(round(sim$get_clock(),2))))
s2 <- plot_state(sim)


sim$update_rates(species = "Clone 1", rates = c(growth = 0.02, death = 0.03))
sim$update_rates(species = "Clone 2", rates = c(growth = 0.4, death = 0.005))
sim$run_up_to_size(species = 'Clone 2', 1000)

t3 <- plot_tissue(sim)+annotate(geom="text",x=900,y=900,label=paste0("Time: ",as.character(round(sim$get_clock(),2))))
s3 <- plot_state(sim)


# Grow 2
sim$update_rates(species = "Clone 1", rates = c(growth = 0.001, death = 0.3))
sim$run_up_to_size(species = 'Clone 2', 2000)
t4 <- plot_tissue(sim)+annotate(geom="text",x=900,y=900,label=paste0("Time: ",as.character(round(sim$get_clock(),2))))
s4 <- plot_state(sim)


# Clone 3
sim$add_mutant(name = "Clone 3", growth_rates = 1, death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.1, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 6000)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.05, death = 0.5))
sim$run_up_to_size("Clone 3", 8000)
t5 <- plot_tissue(sim)+annotate(geom="text",x=900,y=900,label=paste0("Time: ",as.character(round(sim$get_clock(),2))))
s5 <- plot_state(sim)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.01, death = 0.8))
sim$run_up_to_size("Clone 3", 15000) ## orignal 15000
sim$update_rates(species = "Clone 2", rates = c(growth = 0.01, death = 1))
sim$run_up_to_size("Clone 3", 80000)
#t5 <- plot_tissue(sim)+annotate(geom="text",x=900,y=900,label=paste0("Time: ",as.character(round(sim$get_clock(),2))))
#s5 <- plot_state(sim)


# Clone 4
sim$add_mutant(name = "Clone 4", growth_rates = 3 , death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.005, death = 1))
sim$update_rates(species = "Clone 3", rates = c(growth = 0.5, death = 0.02))
sim$mutate_progeny(sim$choose_cell_in("Clone 3"), "Clone 4")
sim$run_up_to_size("Clone 4", 10000) ## original 6000
sim$update_rates(species = "Clone 3", rates = c(growth = 0.02, death = 0.02))
sim$run_up_to_size("Clone 4", 60000) ## orignal 8000
#### try this
#sim$update_rates(species = "Clone 3", rates = c(growth = 0.01, death = 0.05))
#sim$run_up_to_size("Clone 4", 50000)
t6 <- plot_tissue(sim)+annotate(geom="text",x=900,y=900,label=paste0("Time: ",as.character(round(sim$get_clock(),2))))
s6 <- plot_state(sim)

muller <- plot_muller(sim)

print("End simulation")

### Table with info ###
data_final <- inner_join(sim$get_species(),sim$get_counts(),by=c("mutant","epistate")) %>%
                #mutate(percetage=counts/sum(counts)*100) %>%
		select(c("mutant","epistate","growth_rate","death_rate","overall"))
tbl <- ggtexttable(data_final, rows = NULL, theme = ttheme("blank",base_size = 8)) %>%
          tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>%
          tab_add_hline(at.row = 5, row.side = "bottom", linewidth = 3, linetype = 1) %>%
          tab_add_footnote(text = "*Values referring to end of simulation", size = 8, face = "italic")

print("Table written")
### Sampling ###

bboxC <- sim$search_sample(c("Clone 4"= 1000),50, 50)
#bboxC <- sim$search_sample(c("Clone 4" =5000), 100, 100)
#bboxC
#bboxC_lower_corner <- c(180, 300)
#bboxC_upper_corner <- c(280, 400)
#bboxC <- new(TissueRectangle, bboxC_lower_corner, bboxC_upper_corner)
#sim$sample_cells("SPN01_1.3", bboxC_lower_corner, bboxC_upper_corner)
sim$sample_cells("SPN01_1.3", bboxC$lower_corner, bboxC$upper_corner)
print("SPN01_1.3 has been sampled")
t1.2<- plot_tissue(sim)

#bboxA <- sim$search_sample(c("Clone 4" = 1200, 'Clone 3' = 800), 50, 50)
#bboxA <- sim$search_sample(c("Clone 4" = 2500, 'Clone 3' = 2500), 100, 100)
bboxA <- sim$search_sample(c("Clone 4" = 500, 'Clone 3' = 500), 50,50)
bboxA
sim$sample_cells("SPN01_1.1", bboxA$lower_corner, bboxA$upper_corner)
print("SPN01_1.1 has been sampled")
t1.3 <- plot_tissue(sim)

sim
bboxB <- sim$search_sample(c("Clone 3" = 1000),50, 50)
#bboxB <- sim$search_sample(c("Clone 3" = 5000),100, 100)
#bboxB_lower_corner <- c(391,520)
#bboxB_upper_corner <- c(491,620)
#bboxB <- new(TissueRectangle, bboxB_lower_corner, bboxB_upper_corner)
#sim$sample_cells("SPN01_1.2", bboxB_lower_corner, bboxB_upper_corner)
sim$sample_cells("SPN01_1.2", bboxB$lower_corner, bboxB$upper_corner)
t1.1<- plot_tissue(sim)
print("SPN01_1.2 has been sampled")

s = c("Clonal" = "SPN01_1.2", "Clonal" = "SPN01_1.3","Polyclonal" = "SPN01_1.1")
timing = list("T1" = s)

### Use function for plotting different samples
#sampling_plot <- plotting_sample(sim, samples_timing = timing, boxes = box)


# Forest
forest <- sim$get_samples_forest()
forest$save("data/samples_forest_1.sff")

plt_forest <- plot_forest(forest) %>%
	  annotate_forest(forest)+
	   theme(legend.position = "none")
piechart <- plot_state(sim)+theme(legend.position = "none")
timeseries <- plot_timeseries(sim)+
	CNAqc:::my_ggplot_theme()+
	theme(legend.position = "none")

tissue_plots <- list(t1,t2,t3,t4,t5,t6)
state_plots <- list(s1,s2,s3,s4,s5,s6)
plot_tissue_dynamics <- plot_tissue_dynamics(tissue_plots,state_plots,SPN_id="SPN01")
muller <- plot_muller(sim)+
	CNAqc:::my_ggplot_theme()+
	theme(legend.position = "none")
saveRDS(muller,"muller_plot.rds")
page1.1_layout <- "AAAAA\nAAAAA\nAAAAA\nBBBBB\nBBBBB\nCCCDD"
page1.2_layout <- "AAAAA\nBBBBB\nCCCCC\nCCCCC\nCCCCC\nCCCCC"
part1 <-plot_tissue_dynamics # ggplot()
part2 <- ggplot()
part3 <- ggplot()
part4 <- ggplot()
part5 <- ggplot()
### Table with info ###
#data_final <- inner_join(sim$get_species(),sim$get_counts(),by=c("mutant","epistate")) %>%
#	        mutate(percetage=counts/sum(counts)*100)
#tbl <- ggtexttable(data_final, rows = NULL, theme = ttheme("blank",base_size = 8)) %>%
#	  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>%
#	  tab_add_hline(at.row = 5, row.side = "bottom", linewidth = 3, linetype = 1) %>%
#	  tab_add_footnote(text = "*Values referring to end of simulation", size = 8, face = "italic")
page1.1 <- wrap_plots(list(plot_tissue_dynamics,muller,tbl,timeseries),design = page1.1_layout)+
	patchwork::plot_annotation(title= "SPN01 final report",subtitle = "Tissue dynamics")
data_info_samples <- sim$get_samples_info()
tbl_sample <- ggtexttable(data_info_samples, rows = NULL, theme = ttheme("blank",base_size = 8)) %>%
		tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>%
		tab_add_hline(at.row = 4, row.side = "bottom", linewidth = 3, linetype = 1)
ggsave("plot_forest.pdf",plt_forest,width = 250, height = 320, units = "mm", dpi = 300)
#page1.2 <- wrap_plots(list(sampling_plot,tbl_sample,plt_forest),design = page1.2_layout)+
#	        patchwork::plot_annotation(title= "SPN01 final report",subtitle = "Tissue dynamics")
saveRDS(page1.1,"page1.1.rds")
#saveRDS(page1.2,"page1.2.rds")
ggsave("page1.1.pdf",plot=page1.1,width = 250, height = 320, units = "mm", dpi = 300)
#ggsave("page1.2.pdf",plot=page1.2,width = 250, height = 320, units = "mm", dpi = 300)

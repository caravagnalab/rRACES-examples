library(dplyr)
library(hms)
library(ggplot2)
library(patchwork)
library(lubridate)
library(optparse)
source('utils.R')

all_time_sequencing <- tibble()
all_time_merging <- tibble()
all_time_samtools <- tibble()
spns <- paste0('SPN0', 1:5)
for (spn in spns){
  base <- paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/',spn,'/sequencing')
  if (dir.exists(base)){
    print(spn)
    tum <- parse_rds(base, type = 'tumour',merging = FALSE) %>% 
      group_by(coverage, purity, tumour, type, chunk) %>% 
      summarise(elapsed_time_mins = sum(elapsed_time_mins), 
                cpu_time_secs = sum(cpu_time_secs),
                memory_used_MB = sum(memory_used_MB)) %>% 
      mutate(SPN = spn)
    tum_merging <- parse_rds(base, type = 'tumour',merging = TRUE) %>% 
      mutate(SPN = spn)
    
    nor <- parse_rds(base, 'normal', merging = FALSE)  %>% 
      group_by(coverage, purity, tumour, type, chunk) %>% 
      summarise(elapsed_time_mins = sum(elapsed_time_mins), 
                cpu_time_secs = sum(cpu_time_secs),
                memory_used_MB = sum(memory_used_MB)) %>% 
      mutate(SPN = spn)
    nor_merging <- parse_rds(base, type = 'normal',merging = TRUE) %>% 
      mutate(SPN = spn)
    
    time_sequencing <- bind_rows(tum,nor)
    time_merging <- bind_rows(tum_merging,nor_merging)
    all_time_sequencing <- bind_rows(all_time_sequencing, time_sequencing)
    all_time_merging <- bind_rows(all_time_merging, time_merging)
    
    time_t <- parse_out(base, type = 'tumour') %>% 
      mutate(SPN = spn)
    time_n <- parse_out(base, type = 'normal') %>% 
      mutate(SPN = spn)
    time <- bind_rows(time_t, time_n)
    all_time_samtools <- bind_rows(all_time_samtools, time)
  }
}

pl_merging <- all_time_merging %>% 
  ggplot(aes(x = as.factor(type), y = memory_used_MB)) +
  ylab('memory (MB)') + 
  xlab('') +
  geom_boxplot(aes(col = type)) +
  geom_jitter(aes(col = type), height = 0, width = 0.1, size = 0.2) +
  theme_custom() + 
  facet_grid(.~SPN)+
  all_time_merging %>% 
  ggplot(aes(x = as.factor(type), y = hms::as_hms(cpu_time_secs))) +
  ylab('time (H:M:S)') +
  xlab('') + 
  geom_boxplot(aes(col = type)) +
  geom_jitter(aes(col = type), height = 0, width = 0.1, size = 0.2) +
  theme_custom()+ 
  facet_grid(.~SPN) +
  plot_layout(guides = 'collect') + 
  plot_annotation(title = 'Merging')  & theme(legend.position = 'bottom')

pl_sequencing <- all_time_sequencing %>% 
  ggplot(aes(x = as.factor(type), y = memory_used_MB)) +
  ylab('memory (MB)') + 
  xlab('') + 
  geom_boxplot(aes(col = type)) +
  geom_jitter(aes(col = type), height = 0, width = 0.1, size = 0.2) +
  theme_custom() + 
  facet_grid(.~SPN) + 
  all_time_sequencing %>% 
  ggplot(aes(x = as.factor(type), y = hms::as_hms(cpu_time_secs))) +
  ylab('time (H:M:S)') +
  xlab('') + 
  geom_boxplot(aes(col = type)) +
  geom_jitter(aes(col = type), height = 0, width = 0.1, size = 0.2) +
  theme_custom()+ 
  facet_grid(.~SPN) +
  plot_layout(guides = 'collect') + 
  plot_annotation(title = 'Sequencing')  & theme(legend.position = 'bottom')
  
  

pl_samtools <- all_time_samtools %>% 
  ggplot(aes(x = purity, y = memory_MB)) +
  geom_boxplot(aes(col = type)) +
  geom_jitter(aes(col = type),height = 0, width = 0.1, size = 0.5) +
  ylab('memory (MB)') +
  xlab('') + 
  facet_grid(step~SPN, scales = 'free_y') + 
  theme_custom()  + 

  all_time_samtools %>% 
  ggplot(aes(x = purity, y = as.POSIXct(usr_time))) +
  geom_boxplot(aes(col = type)) +
  scale_y_datetime(
    date_labels = "%H:%M"
  ) + 
  geom_jitter(aes(col = type), height = 0, width = 0.1, size = 0.5) +
  facet_grid(step~SPN, scales = 'free_y') + 
  ylab('usr_time (H:M)') +
  xlab('') + 
  theme_custom() +
  plot_layout(nrow = 2, guides = 'collect') + 
  plot_annotation(title = 'samtools') & theme(legend.position = 'bottom')

patchwork::wrap_plots(pl_sequencing, pl_merging,pl_samtools, nrow = 3, design = 'A\nB\nC\nC\nC')+plot_layout(guides = "collect")
#ggsave(filename ='plot_benchmark.png', dpi = 300, units = 'in', width = 14, height = 16)
#ggsave(filename = paste0(param$outputdir, '/plot_benchmark.png'), dpi = 300, units = 'in', width = 10, height = 12)




library(dplyr)
library(hms)
library(ggplot2)
library(patchwork)
library(lubridate)
library(optparse)

option_list <- list(
  make_option(c("-i", "--inputdir"), type="character", default='',
              help="Input directory"),
  make_option(c("-o", "--outputdir"), type="character", default='',
              help="Output directory")
)
param <- parse_args(OptionParser(option_list=option_list))

base <- param$inputdir

parse_rds <- function(dir, type){
  df <- tibble()
  
  if (type == 'normal'){
    path <-  list(paste(dir, type, 'purity_1/data/resources/', sep = '/'))
  } else if (type == 'tumour'){
    path <- list.dirs(paste(dir, type,sep = '/'), recursive = F)
    path <- paste(path, 'data/resources/', sep = '/')
  }
  
  for (p in path){
    purity <- strsplit(p, split = '/') %>% unlist()
    purity <- purity[[11]]
    
    files <- list.files(p)
    for (f in files){
      chunk <- strsplit(gsub(pattern = '_', replacement = '.', x = f), '\\.')  %>% unlist()
      chunk <- chunk[length(chunk)-1]
      
      data <- readRDS(paste0(p, f)) %>% tibble()
      data <- data %>% mutate(chunk = chunk) %>% mutate(purity = purity)
      df <- bind_rows(df, data)
    }
  }
  df <- df %>% mutate(type = type)
  return(df)
}
tum <- parse_rds(base, type = 'tumour') %>% 
  group_by(coverage, purity, tumour, type, chunk) %>% 
  summarise(elapsed_time_mins = sum(elapsed_time_mins), 
            cpu_time_secs = sum(cpu_time_secs),
            memory_used_MB = sum(memory_used_MB))

nor <- parse_rds(base, 'normal')  %>% 
  group_by(coverage, purity, tumour, type, chunk) %>% 
  summarise(elapsed_time_mins = sum(elapsed_time_mins), 
            cpu_time_secs = sum(cpu_time_secs),
            memory_used_MB = sum(memory_used_MB))
time_sequencing <- bind_rows(tum,nor)

pl_sequencing <- time_sequencing %>% 
  ggplot(aes(x = as.factor(purity), y = memory_used_MB)) +
  ylab('memory (MB)') + 
  xlab('purity') + 
  geom_violin() +
  geom_jitter(aes(col = type), height = 0, width = 0.1, size = 0.5) +
  theme_custom() + 
time_sequencing %>% 
  ggplot(aes(x = as.factor(purity), y = hms::as_hms(cpu_time_secs))) +
  ylab('time (H:M:S)') +
  xlab('purity') + 
  geom_violin() +
  geom_jitter(aes(col =type), height = 0, width = 0.1, size = 0.5) +
  theme_custom()+ 
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  #plot_annotation(title = c('simulate_seq')) 
  
  


parse_out <- function(dir, type){
  if (type == 'normal'){
    path <-  list(paste(dir, type, 'purity_1/TIME/', sep = '/'))
  } else if (type == 'tumour'){
    path <- list.dirs(paste(dir, type,sep = '/'), recursive = F)
    path <- paste(path, 'TIME/', sep = '/')
  }
  
  df <- tibble()
  
  for (p in path){
    purity <- strsplit(p, split = '/') %>% unlist()
    purity <- purity[[11]] %>% strsplit('_') %>% unlist() %>% tail(n = 1)
    
    all_steps <- c('samtools_merge', 'fastq', 'samtools_split')
    for (step in all_steps){
        step_df <- tibble()
        
        all_files <- list.files(path = p) 
        files <- grep(step, all_files, value = TRUE)
        
        for (f in files){
          chunk <- strsplit(f, '_') %>% unlist()
          
          if (step == 'fastq'){
            chunk <- chunk[length(chunk)-2]
          } else if (step == 'samtools_merge'){
            chunk <- chunk[length(chunk)-2]
          } else{
            chunk <- strsplit(chunk[length(chunk)], '\\.') %>% unlist() %>% head(1)
          }
          
          data <- readLines(paste0(p, f))
          usr_time_sec <- strsplit(data[[2]] , ':') %>% unlist() %>% tail(1) %>% trimws() %>% as.numeric()
          usr_time <- strsplit(data[[2]] , ':') %>% unlist() %>% tail(1) %>% trimws() %>% as.numeric() %>% as_hms()
          sys_time_sec <- strsplit(data[[3]] , ':') %>% unlist() %>% tail(1) %>% trimws() %>% as.numeric()
          sys_time <- strsplit(data[[3]] , ':') %>% unlist() %>% tail(1) %>% trimws() %>% as.numeric() %>% as_hms()
          elapsed_time <- strsplit(gsub(pattern = '\\(h:mm:ss or m:ss\\):', replacement = '_', x = data[[5]]), '_') %>% unlist() %>% tail(1) %>% trimws()
          mem_kb <- strsplit(data[[10]] , ':') %>% unlist() %>% tail(1) %>% trimws()
          
          tmp_df <- tibble(usr_time = usr_time,
                           usr_time_sec = usr_time_sec,
                           sys_time = sys_time,
                           sys_time_sec = sys_time_sec,
                           elapsed_time = elapsed_time,
                           memory_GB = as.numeric(mem_kb)/(1024^2),
                           memory_MB = as.numeric(mem_kb)/1024,
                           chunk = chunk,
                           purity = purity,
                           type = type)
          step_df <- bind_rows(step_df, tmp_df)
        }
        step_df <- step_df %>% mutate(step = step) 
        df <- bind_rows(df, step_df)
      }
  }
  return(df)
}

time_t <- parse_out(base, type = 'tumour')
time_n <- parse_out(base, type = 'normal')
time <- bind_rows(time_t, time_n)

pl_samtools <- time %>% 
  ggplot(aes(x = purity, y = memory_MB)) +
  geom_violin() +
  geom_jitter(aes(col = type),height = 0, width = 0.1, size = 0.5) +
  ylab('memory (MB)') +
  facet_wrap(.~step, scales = 'free_y') + 
  theme_custom()  + 

time %>% 
  ggplot(aes(x = purity, y = as.POSIXct(usr_time))) +
  geom_violin() +
  scale_y_datetime(
    date_labels = "%H:%M"
  ) + 
  geom_jitter(aes(col = type), height = 0, width = 0.1, size = 0.5) +
  facet_wrap(.~step, scales = 'free_y') + 
  ylab('usr_time (H:M)') +
  theme_custom() +
  plot_layout(nrow = 2) + 
  plot_annotation(title = 'samtools')

patchwork::wrap_plots(pl_sequencing, pl_samtools, nrow = 2, design = 'A\nB\nB\nB')
ggsave(filename = paste0(param$outputdir, '/plot_benchmark.png'), dpi = 300, units = 'in', width = 10, height = 12)


# plot total time
# sm_time_all <- time %>% group_by(chunk, step) %>% summarise(max_time = max(usr_time)) %>% group_by(chunk) %>% summarize(tot_time = as_hms(sum(max_time)))
# seq_time_t <- tum %>% select(time, chunk, type)
# 
# full_join(sm_time_all, seq_time_t, by = join_by(chunk)) %>% 
#   mutate(final_time = as_hms(tot_time+ time)) %>% 
#   ggplot(aes(x = type, y = as.POSIXct(final_time))) +
#   geom_violin() +
#   scale_y_datetime(
#     date_labels = "%H:%M"
#   ) + 
#   geom_jitter(height = 0, width = 0.1, size = 0.5) +
#   ylab('total time (H:M)') +
#   theme_custom() 






theme_custom <- function() {
  theme_bw(base_size = 13) + # Start with a minimal theme
    theme(
      # Background color
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Title and axis labels
      plot.title = element_text(size = 15, color = "#333333"), #hjust = 0.5),
      plot.subtitle = element_text(size = 13, color = "#444444"),
      axis.title = element_text(size = 15, color = "#444444"),
      axis.text = element_text(size = 14, color = "#555555"),
      
      # Grid lines
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_line(color = "grey90", size = 0.5),
      
      # Legend
      legend.background = element_rect(fill = "white", color = "NA"),
      legend.text = element_text(size = 14, color = "#444444"),
      legend.title = element_text(size = 14, color = "#444444"),
      legend.position = 'bottom',
      legend.box.spacing = unit(-0.5, "pt"),
      
      # Margins and padding
      plot.margin = margin(10, 10, 10, 10),
      
      # Optional: Customize axis lines or ticks
      axis.line = element_line(color = "grey70"),
      axis.ticks = element_line(color = "grey70")
    )
}
                                       

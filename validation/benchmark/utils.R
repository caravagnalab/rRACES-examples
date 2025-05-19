parse_rds <- function(dir, type, merging){
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
    if (merging){
      files <- list.files(p,pattern = "seq_results_merging")
      for (f in files){
        coverage <- strsplit(gsub(pattern = '_', replacement = '.', x = f), '\\.')  %>% unlist()
        coverage <- coverage[length(coverage)-1]
        data <- readRDS(paste0(p, f)) %>% tibble()
        data <- data %>% mutate(coverage = coverage) %>% mutate(purity = purity)
        df <- bind_rows(df, data)
      }
    } else{
      files <- list.files(p,pattern = "seq_results_resources")
      for (f in files){
        chunk <- strsplit(gsub(pattern = '_', replacement = '.', x = f), '\\.')  %>% unlist()
        chunk <- chunk[length(chunk)-1]
        
        data <- readRDS(paste0(p, f)) %>% tibble()
        data <- data %>% mutate(chunk = chunk) %>% mutate(purity = purity)
        df <- bind_rows(df, data)
      }
    }     
  }
  df <- df %>% mutate(type = type)
  return(df)
}




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
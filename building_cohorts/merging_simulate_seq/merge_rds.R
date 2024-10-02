library(dplyr)
library(rRACES)
library(optparse)
option_list <- list(
  make_option(c("--lot_range"), type = "character", default = NULL,
              help = "Range in which to split the different lots", metavar = "character"),
  make_option(c("--lot_type"), type = "character", default = NULL,
              help = "Type of lot, either single or final", metavar = "character"),
  make_option(c("--rds_dir"), type = "character", default = NULL,
              help = "Path to the directory where simulate_seq.rds are stored", 
              metavar = "character"),
  make_option(c("--chrom"), type = "numeric", default = NULL,
              help = "Chromosome number", 
              metavar = "numeric")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

lot_range <- opt$lot_range
lot <- opt$lot_type
dir<-opt$rds_dir
chrom<-opt$chrom

# lot_range <- args[1]
# dir <- args[3]
# dir <- "10X_0.9p/data"
splited_range <- strsplit(lot_range, ":", fixed = TRUE)[[1]]
lot_start <- as.numeric(splited_range[1]) ## start of the interval
lot_end <- as.numeric(splited_range[2]) ## end of the interval
# lot <- args[2]
setwd(dir)
dir.create(path = chrom)
data <- list()
if (lot == "single") {
	  rds_files <- list.files(path = ".", pattern = "seq_results_SPN01_")
  for (i in c(lot_start:lot_end)){
    print(rds_files[i])
    d <- readRDS(rds_files[i]) %>% 
      dplyr::filter(chr==chrom)
    data[[i]] <- seq_to_long(d) %>%
      dplyr::select(!VAF)
		print(i)
	}
} else {
	  rds_files <- list.files(path = ".", pattern = "seq_results_lot_")
  for (i in c(lot_start:lot_end)){
	      print(rds_files[i])
	      data[[i]] <- readRDS(rds_files[i]) %>% 
	        dplyr::filter(chr==chrom)
    }
}
print("Combining dataframes ...")
combined_df <- bind_rows(data)


print("Summing up NV and DP...")
# Group by 'mutID' and sum 'NV' and 'DP'
# result <- combined_df %>%
# 	  group_by(chr, from, sample_name, ref, alt, causes, classes) %>%
# 	    summarise(
# 		          NV = sum(NV, na.rm = TRUE),
# 			      DP = sum(DP, na.rm = TRUE)
# 			    )

combined_df_long_DP <-combined_df %>% 
  group_by(chr, from, sample_name,classes,causes) %>%
  summarise(
    DP = sum(DP, na.rm = TRUE))

combined_df_long_NV <-combined_df %>% 
  group_by(chr, from, ref,alt,sample_name,classes,causes) %>%
  summarise(
    NV = sum(NV, na.rm = TRUE))

result <- full_join(x = combined_df_long_DP,
                                  y = combined_df_long_NV,
                                  by = c("chr","from","sample_name","classes","causes")) %>% 
  mutate("VAF"=NV/DP)



if (lot == "final") {
  saveRDS(result, file = paste0(chrom,"/","seq_results_final_",chrom,".rds"))
} else {
  saveRDS(result, file = paste0(chrom,"/","seq_results_lot_", chrom,"_",lot_range, ".rds"))
}

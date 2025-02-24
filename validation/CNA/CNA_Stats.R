# Load necessary libraries
library(dplyr)
library(GenomicRanges)
library(rRACES)
library(optparse)
library(tidyr)

option_list <- list(
  make_option(c("--phylo_forest_path"), type = "character", default = NULL,
              help = "Path to the phylogenetic forest file", metavar = "character"),
  make_option(c("--sample_id"), type = "character", default = NULL,
              help = "Name of the rRACES sample", metavar = "character"),
  make_option(c("--segment_file"), type = "character", default = NULL,
              help = "Path to the segment file", metavar = "character"),
  make_option(c("--overlap_threshold"), type = "numeric", default = NULL,
              help = "Threshold for overlapping", metavar = "numeric")
#  make_option(c("--samples"), type = "character", default = NULL,
#              help = "Comma-separated list of samples", metavar = "character"),
#  make_option(c("--samples_vc"), type = "character", default = NULL,
#              help = "Comma-separated list of samples for variant calling", metavar = "character")

)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

phylo_forest_path <- opt$phylo_forest_path
sample_id <- opt$sample_id
segment_file<-opt$segment_file
overlap_threshold<-opt$overlap_threshold

phylo_forest <- load_phylogenetic_forest(phylo_forest_path)
bulk_sample <- phylo_forest$get_bulk_allelic_fragmentation(sample_id)
ascat_seg <- read.table(segment_file,header=T)



# Function to calculate overlap ratio between two segments
calculate_overlap_ratio <- function(start1, end1, start2, end2) {
	overlap_start <- max(start1, start2)
  	overlap_end <- min(end1, end2)
  	overlap_length <- max(0, overlap_end - overlap_start + 1) # overlap length

    	union_length <- max(end1, end2) - min(start1, start2) + 1 # union length

      	return(overlap_length / union_length) # overlap ratio
}

# Add segment IDs to both ascat_seg and bulk_sample
ascat_seg <- ascat_seg %>%
	mutate(segment_id = paste0("ascat_", chr,":",startpos,":",endpos,"_",nMajor,":",nMinor))

bulk_sample <- bulk_sample %>%
	mutate(segment_id = paste0("bulk_", chr,":",begin,":",end,"_",major,":",minor))

# Initialize a dataframe to store pairwise comparison results
pairwise_comparison <- data.frame(
	ascat_segment_id = character(),
	bulk_segment_id = character(),
	overlap_ratio = numeric(),
	nMajor_match = logical(),
	nMinor_match = logical(),
	stringsAsFactors = FALSE
)

# Loop through all chromosomes
for (chromosome in unique(ascat_seg$chr)) {
	# Filter segments by chromosome
	ascat_chr <- ascat_seg %>% filter(chr == chromosome)
      	bulk_chr <- bulk_sample %>% filter(chr == chromosome)
        
        # Compare all possible pairs between ascat_seg and bulk_sample for this chromosome
        for (i in 1:nrow(ascat_chr)) {
		for (j in 1:nrow(bulk_chr)) {
			# Calculate overlap ratio
			overlap_ratio <- calculate_overlap_ratio(ascat_chr$startpos[i], ascat_chr$endpos[i],
							         bulk_chr$begin[j], bulk_chr$end[j])
            
        # Check if nMajor and nMinor match
        nMajor_match <- ascat_chr$nMajor[i] == bulk_chr$major[j]
	nMinor_match <- ascat_chr$nMinor[i] == bulk_chr$minor[j]
	          
	# Store the comparison result
	pairwise_comparison <- rbind(pairwise_comparison, 
			  	data.frame(ascat_segment_id = ascat_chr$segment_id[i],
					bulk_segment_id = bulk_chr$segment_id[j],overlap_ratio = overlap_ratio,
					nMajor_match = nMajor_match,nMinor_match = nMinor_match,															                                               		      stringsAsFactors = FALSE))
		}
        }
    }
pairwise_comparison <- pairwise_comparison %>%
	separate(ascat_segment_id, into = c("caller", "ascat_segment", "ascat_karyotype"), sep = "_") %>%
	separate(bulk_segment_id, into = c("races", "races_segment", "races_karyotype"), sep = "_")


# Set an overlap threshold (e.g., 0.5)
#overlap_threshold <- 0.5

# Determine the best match for each segment in ascat_seg
best_matches <- pairwise_comparison %>%
		  filter(overlap_ratio >= overlap_threshold) %>%
		  group_by(ascat_segment) %>%
		  top_n(1, wt = overlap_ratio)
saveRDS(best_matches,paste0(sample_id,"_best_matches.rds"))
# Count correct matches based on copy number
correct_segments <- sum(best_matches$nMajor_match & best_matches$nMinor_match)

# Calculate precision, recall, and F1-score
total_predicted_segments <- nrow(ascat_seg)
total_true_segments <- nrow(bulk_sample)
precision <- correct_segments / total_predicted_segments
recall <- correct_segments / total_true_segments
f1_score <- 2 * ((precision * recall) / (precision + recall))
cat("Correct Segments:", correct_segments, "\n")
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1-score:", f1_score, "\n")

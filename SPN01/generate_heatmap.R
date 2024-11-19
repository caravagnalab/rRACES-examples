library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)

########################################
######## Some utils functions ##########
########################################

get_segment_info = function(data, chr, sample, new_from, new_to, keep_columns){
  data %>%
    dplyr::filter(sample_id == sample,
                  chr == !!chr,
                  from <= new_from,
                  to >= new_to) %>%
    select(all_of(keep_columns))
}

average_in_window = function(window, gr, v, method = "weighted", empty_v = NA) {
  
  if(missing(v)) v = rep(1, length(gr))
  if(is.null(v)) v = rep(1, length(gr))
  if(is.atomic(v) && is.vector(v)) v = cbind(v)
  
  v = as.matrix(v)
  if(is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }
  
  if(length(empty_v) == 1) {
    empty_v = rep(empty_v, ncol(v))
  }
  
  u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))
  
  mtch = as.matrix(findOverlaps(window, gr))
  intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
  w = width(intersect)
  v = v[mtch[,2], , drop = FALSE]
  n = nrow(v)
  
  ind_list = split(seq_len(n), mtch[, 1])
  window_index = as.numeric(names(ind_list))
  window_w = width(window)
  
  if(is.character(v)) {
    for(i in seq_along(ind_list)) {
      ind = ind_list[[i]]
      if(is.function(method)) {
        u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
      } else {
        tb = tapply(w[ind], v[ind], sum)
        u[window_index[i], ] = names(tb[which.max(tb)])
      }
    }
  } else {
    if(method == "w0") {
      gr2 = reduce(gr, min.gapwidth = 0)
      mtch2 = as.matrix(findOverlaps(window, gr2))
      intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])
      
      width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
      ind = unique(mtch2[, 1])
      width_setdiff = width(window[ind]) - width_intersect
      
      w2 = width(window[ind])
      
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
        u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
      }
      
    } else if(method == "absolute") {
      for(i in seq_along(ind_list)) {
        u[window_index[i], ] = colMeans(v[ind_list[[i]], , drop = FALSE])
      }
      
    } else if(method == "weighted") {
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        or <- colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
        u[window_index[i], ] = round(or,0)
      }
    } else {
      if(is.function(method)) {
        for(i in seq_along(ind_list)) {
          ind = ind_list[[i]]
          u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
        }
      } else {
        stop("wrong method.")
      }
    }
  }
  
  return(u)
}


########################################
bulk_cell_cna <- readRDS("bulk_cell_cna.rds")
x <- bulk_cell_cna %>% 
  dplyr::mutate(chr=paste0("chr",chr)) %>% 
  dplyr::rename("from"="begin") %>% 
  dplyr::rename("to"="end") %>% 
  dplyr::rename(Major=major) %>% 
  dplyr::rename(sample_id=sample) %>% 
  mutate(segment_id=paste(chr,from,to,sep=":"))


out = lapply(x$chr %>% unique(), function(chr) {
  
  cli::cli_alert_info("Iterating on {.field {chr}}")
  cli::cli_h2("Iterating on {.field {chr}}")
  
  old_segments = sapply(x$sample_id %>% unique(), function(s) {
    x %>%
      dplyr::filter(sample_id == s) %>%
      dplyr::filter(chr == !!chr) %>%
      dplyr::pull(segment_id) %>%
      unique() %>%
      length()
  })
  
  cli::cli_alert_info("Number of original segments in individual {.cls CNAqc} objects in {.field {chr}}:")
  cli::cli_ul(paste(names(old_segments), old_segments, sep = " = "))
  cat("\n")
  
  # Chromosome-specific new breakpoints
  new_breakpoints = c(
    x %>%
      dplyr::filter(chr == !!chr) %>%
      dplyr::pull(from),
    x %>%
      dplyr::filter(chr == !!chr) %>%
      dplyr::pull(to)) %>% 
    unique() %>%
    sort()
  
  cli::cli_alert_info("Found {.field {length(new_breakpoints)}} breakpoints")
  #cat("\n")
  
  #  Separate new breakpoints into segment from and to values
  new_from = new_breakpoints[ !new_breakpoints == dplyr::last(new_breakpoints)] # last element will not be included in the from column 
  new_to = new_breakpoints[!new_breakpoints == dplyr::first(new_breakpoints)] # first element will not be included in the to column 
  
  # iterate on the new breakpoints to subset the cna piled up 
  lapply(new_breakpoints[-1] %>% seq_along(), function(i) {
    
    # iterate for each sample 
    lapply(x$sample_id %>% unique(), function(s) {
      
      # define which columns must be kept in the new table
      not_wanted = c("from", "to", "length", "size", "segment_id", "chr", "sample_id", "n")
      wanted = setdiff(colnames(x), not_wanted)
      
      # get the information for the sample in the new segment 
      tmp = get_segment_info(x, 
                             chr = chr, 
                             sample = s, 
                             new_from = new_from[i], 
                             new_to = new_to[i], 
                             keep_columns = wanted)
      
      # do some checking on the result
      if (nrow(tmp) == 0) { # there is no information on copy number on the new segment: insert NA as value of all the columns, except from, to, segment_id and sample_id
        
        tmp_v2 = rep(NA, ncol(tmp))  
        names(tmp_v2) = colnames(tmp)
        
        tmp = tmp_v2 %>% 
          tibble::as_tibble_row()
      }
      
      # create a tibble with the information on the new breakpoints and include the previously retrieved information 
      tidyr::tibble(chr = chr, 
                    from = new_from[i], 
                    to = new_to[i], 
                    sample_id = s) %>% 
        dplyr::bind_cols(tmp) %>% 
        dplyr::mutate(segment_id = paste(chr, from, to, sep = ":")) 
      
    }) %>% do.call(bind_rows, .)
  }) %>% do.call(bind_rows, .)
})  %>% do.call(bind_rows, .) # create a unique big tibble with the new segmentation



wide_df <- out %>%
  na.omit() %>%
  select(chr, from, to, cell_id, total_cn) %>% # Select only relevant columns
  pivot_wider(names_from = cell_id, values_from = total_cn,values_fill = list(total_cn = 0),values_fn = list(total_cn = sum)) %>% as.data.frame() # Reshape to wide format

wide_df_gr1 <- GRanges(seqnames = wide_df[, 1], ranges = IRanges(wide_df[, 2], wide_df[, 3]))



chr_df = read.chromInfo()$df
chr_df = chr_df[chr_df$chr %in% paste0("chr", 1:22), ]
chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
chr_gr

library(EnrichedHeatmap)
chr_window = makeWindows(chr_gr, w = 1e6)
chr_window

chr = as.vector(seqnames(chr_window))
chr_level = paste0("chr", 1:22)
chr = factor(chr, levels = chr_level)


num_mat_mine = average_in_window(chr_window, wide_df_gr1, wide_df[, -(1:3)])
num_mat_mine[is.na(num_mat_mine)] <- 0

#subgroup = 
## c(unlist(out[,"mutant"]))
#subgroup = out %>% select(cell_id,mutant) %>% na.omit() %>% unique()
subgroup = out %>% select(cell_id,mutant,sample_id) %>% na.omit() %>% unique()


subgroup <- subgroup %>%
  group_by(sample_id) %>%
  mutate(
    clonal_status = ifelse(n_distinct(mutant) > 1, "polyclonal", "monoclonal")
  ) %>%
  ungroup()

######## for karotype
wide_df_kar <- out %>%
  na.omit() %>% 
  mutate(cell_id=as.character(cell_id)) %>% 
  select(chr, from, to, cell_id, karyotype) %>% # Select only relevant columns
  pivot_wider(names_from = cell_id, values_from = karyotype) %>% as.data.frame() # Reshape to wide format
char_mat = NULL
n_cells <- ncol(wide_df_kar)-3
for(i in 1:n_cells) {
  reg = wide_df_kar[,1:3]
  j=i+3
  reg = cbind(reg,wide_df_kar[,j])
  # bed = bed[sample(nrow(bed), 20), , drop = FALSE]
  gr_cnv = GRanges(seqnames = reg[, 1], ranges = IRanges(reg[, 2], reg[, 3]))
  
  char_mat = cbind(char_mat, average_in_window(chr_window, gr_cnv, reg[, 4]))
  print(i)
}


colors = c("0" = "steelblue", "1" = "skyblue1","2" = "honeydew3","3"="peachpuff","4"="tan1",
           "5"="tan4","6"="red4","7"="violetred4","8"="violetred1")


col_clone = RColorBrewer::brewer.pal(4, "Paired")[3:4]
names(col_clone)=unique(subgroup$mutant)
col_sample <- RColorBrewer::brewer.pal(3, "Dark2")
names(col_sample) = unique(subgroup$sample_id)
col_clonality = c("monoclonal" = "#EEDD82","polyclonal"="lightgoldenrod4")
row_ha = rowAnnotation(clone = subgroup$mutant,
                           sample = subgroup$sample_id,
                           clonality=subgroup$clonal_status,
                           ploidy=(rowMeans(t(num_mat_mine))),
                       col=list(clone=col_clone,sample=col_sample,clonality=col_clonality))



colors <- rep(c("gray", "black"), length.out = length(levels(chr)))

# Create a named vector to assign colors to each chromosome level
color_map <- setNames(colors, levels(chr))

# Apply the color map to the factor to get the color for each chromosome in `chr`
chr_colors <- color_map[chr]
column_ha = HeatmapAnnotation(chromosome = chr, col=list(chromosome=chr_colors),show_legend = F)





pdf("heatmap.pdf",width=20,height=20)
Heatmap(t(char_mat), name = "CNV", col = c("1:1" = "#228B22CC", "1:0" = "steelblue",
                                           "0:0"="darkblue","2:0"="turquoise4",
                                           "2:1"="#FFA500CC","2:2"="firebrick3","Other"="grey"),
        top_annotation = column_ha,
        right_annotation = row_ha,
        column_split = chr, cluster_columns = FALSE, show_row_dend = FALSE,column_title_gp = gpar(fontsize = 10), border = TRUE,
        column_gap = unit(0, "points"),
        column_title = ifelse(1:22 %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n")),
        heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"))

dev.off()






# library(ComplexHeatmap)
# pdf("heatmap.pdf",width=20,height=20)
# ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
# ht_list =Heatmap(t(num_mat_mine), name = "mat",col=colors,
#                  column_split = chr, cluster_columns = FALSE, show_row_dend = FALSE,
#                  row_split = subgroup$mutant, cluster_row_slices = FALSE,
#                  row_title = "numeric matrix",
#                  left_annotation = rowAnnotation(clone = subgroup$mutant, sample=subgroup$sample_id,show_annotation_name = FALSE,
#                                                  annotation_legend_param = list(
#                                                    subgroup = list(direction = "horizontal", title_position = "lefttop", nrow = 1))),
#                  column_title_gp = gpar(fontsize = 10), border = TRUE,
#                  column_gap = unit(0, "points"),
#                  column_title = ifelse(1:22 %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n")),
#                  heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"))+
#   Heatmap(t(char_mat), name = "CNV", col = c("1:1" = "#228B22CC", "1:0" = "steelblue",
#                                              "0:0"="darkblue","2:0"="turquoise4",
#                                              "2:1"="#FFA500CC","2:2"="firebrick3","Other"="grey"),
#           border = TRUE, column_title = "character matrix",
#           column_split = chr, cluster_columns = FALSE, show_row_dend = FALSE)
# draw(ht_list, merge_legend = TRUE)
# dev.off()

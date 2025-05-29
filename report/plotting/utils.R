library(ProCESS)
library(dplyr)
library(patchwork)
library(ggplot2)
#library(CNAqc)
#source("/orfeo/cephfs/scratch/cdslab/ggandolfi/prj_races/ProCESS-examples/SPN01/scripts/my_functions/plot_genome_wide.R")
#'
#' Convert realtive coordinates to absolute coordinates
#'
#' This function takes CNAqc bject and converts relative coordinates
#' to absolute coordinates given an input reference genome.
#'
#' @param x CNAqc object
#' @param ref name of the reference genome
#' @return x CNAqc object with absoluute coordinates
#'
#' @export


relative_to_absolute_coords_pos = function(x, ref = "GRCh38") {
  reference_genome = CNAqc:::get_reference(ref)
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  x = x %>%
      dplyr::rename(pos=from) %>%
      dplyr::mutate(chr = paste0("chr", chr)) %>%
      dplyr::mutate(pos = pos + vfrom[chr])
  return(x)
}


#'
#' Generates CNAqc object given simulate_seq() dataframe and
#' phylogenetic forest object. The final result will be a CNAqc 
#' object.
#'

races2cnaqc <- function(seq_results,phylo_forest,sample_id,ref,purity){
  ref_path <- phylo_forest$get_reference_path()
  driver_table_path <- gsub(pattern = "reference.fasta",replacement = "drivers.txt",x = ref_path)
  driver_table <-  read.csv(driver_table_path,header=T,sep="\t")
  known_drivers <- driver_table %>%
    dplyr::mutate(chr=as.character(chr)) %>%
    dplyr::rename(driver_label=driver_code)
  bulk <- phylo_forest$get_bulk_allelic_fragmentation(sample_id)
  cna <- bulk %>% dplyr::rename("Major"=major,"from"=begin,"to"=end) %>%
    dplyr::filter(ratio>=0.08)
  cna <- cna %>%
    group_by(chr, from, to) %>%            # Group by chr, begin, and end
    arrange(desc(ratio)) %>%                 # Arrange by 'ratio' in descending order
    mutate(ccf_label = paste0("ccf_", rank(-ratio)))  # Assign rank based on 'ratio'

  mutations <- ProCESS::seq_to_long(seq_results) %>%
    dplyr::filter(sample_name==sample_id & classes!="germinal") %>%
    dplyr::filter(VAF!=0) %>% mutate(is_driver=FALSE) %>%
    left_join(known_drivers,by=c("chr","from","to","ref","alt")) %>%
    dplyr::mutate(
      is_driver = ifelse(!is.na(driver_label), TRUE,
                         ifelse(is.na(driver_label) & classes == "driver", TRUE, FALSE)),
      driver_label = ifelse(is.na(driver_label) & classes == "driver",
                            paste(chr, from, ref, alt, sep = ":"), driver_label))
  x <- CNAqc::init(mutations = mutations,cna = cna,
                   purity = purity,sample = sample_id,
                   ref = ref)
  return(x)
}


get_classes_colors <- function(classes){
  color = c(`driver` = "firebrick4",`passenger` = ggplot2::alpha("tan2",0.4),`pre-neoplastic` = "cornflowerblue",
            `germinal` = "darkolivegreen")
  missing = setdiff(names(color), classes)
  nmissing = length(missing)
  c(color, CNAqc:::nmfy(missing, rep("gray", nmissing)))
}

# get_clone_palette <- function(sample_forest){
#   clones <- sample_forest
# }


get_legend <- function(col_palette){
  df <- data.frame(type = names(col_palette), color = col_palette)
  p <- ggplot(df, aes(x = type, fill = type)) +
    geom_bar() +
    scale_fill_manual(values = col_palette) +
    theme_void() +  # Remove axes and background
    guides(fill = guide_legend(title = "Classes & Causes"))

  legend_plot <- ggpubr::get_legend(p,position = "right")
  pl <- ggpubr::as_ggplot(legend_plot)
  return(pl)
}


squareplot = function(seq_res, samples_list,chrom)
{
  row_plots = NULL
  for (s in seq(samples_list))
  {
    sn = samples_list[s]
    s_seq <- seq_res %>% filter(classes!="germinal")
    s_seq_long <- s_seq %>% ProCESS::seq_to_long()
    plot_vaf <- s_seq_long %>%
      filter(sample_name==sn & chr==chrom) %>%
      filter(VAF!=0) %>%
      ggplot(aes(x=VAF)) +geom_histogram(binwidth = 0.01) +
      xlim(c(0,1))+
      ggplot2::ggtitle(label = sn) +
      CNAqc:::my_ggplot_theme()



    mb = list(plot_vaf+ labs(title = sn) )

    idx_pre = 1:s
    idx_post = s:length(samples_list)

    pl_r = pl_l = NULL
     
    #palette <- RColorBrewer::brewer.pal(n = length(unique(s_seq$causes)), name = "Set3")
    #col_causes <- setNames(palette, unique(s_seq$causes))
    #cols_causes <- ProCESS:::get_colors_for(unique(s_seq$causes))
    #col_classes <- c("passenger" = "#CCCC99",
    #                 "pre-neoplastic" = "#006699",
    #                 "driver" = "#990033")
    #cols <- c(col_causes,col_classes)

    if (length(idx_pre) > 1)
      pl_r = lapply(setdiff(idx_pre, s), function(x) {
        s_sn <- s_seq_long %>% filter(sample_name==sn & chr==chrom)
        s_sn_x <- s_seq_long %>% filter(sample_name==samples_list[x] & chr==chrom)
        joined <- full_join(s_sn,s_sn_x,by=c("chr","from","ref","alt","to","causes","classes"))
        plot <- joined %>% ggplot(aes(x=VAF.x,y=VAF.y,col=classes)) + geom_point() +
          CNAqc:::my_ggplot_theme()
          #scale_color_manual(values = col_classes)
        plot + ggplot2::geom_point(alpha = 0.7) +
          ggplot2::xlim(c(-0.01, 1.01)) +
          ggplot2::ylim(c(-0.01, 1.01)) +
          ggplot2::labs(x = sn, y = samples_list[x])+
          ggplot2::theme(legend.position = "none")
      })

    if (length(idx_post) > 1)
      pl_l = lapply(setdiff(idx_post, s), function(x) {
        s_sn <- s_seq_long %>% filter(sample_name==sn & chr==chrom)
        s_sn_x <- s_seq_long %>% filter(sample_name==samples_list[x] & chr==chrom)
        joined <- full_join(s_sn,s_sn_x,by=c("chr","from","ref","alt","to","causes","classes"))
        plot <- joined %>% ggplot(aes(x=VAF.x,y=VAF.y,col=causes)) + geom_point() +
          CNAqc:::my_ggplot_theme()
          #scale_color_manual(values = col_causes)
        plot + ggplot2::geom_point(alpha = 0.7) +
          ggplot2::xlim(c(-0.01, 1.01)) +
          ggplot2::ylim(c(-0.01, 1.01)) +
          ggplot2::labs(x = sn, y = samples_list[x])+
          ggplot2::theme(legend.position = "none")
      })

    plotlist = append(append(pl_r, mb), pl_l)
    row_plot = patchwork::wrap_plots(plotlist)+
      patchwork::plot_layout(guides = "collect",ncol = length(pl_r) + length(pl_l) + 1,nrow = 1)
    row_plots = append(row_plots, list(row_plot))
  }
  #pl <- get_legend(cols)
  patchwork::wrap_plots(row_plots)+
    patchwork::plot_layout(design = "AAAA\nBBBB\nCCCC")+
    patchwork::plot_annotation(title = paste0("Chromosome ", chrom))
}

get_clone_map <- function(sample_forest){
  n_clones <- nrow(sample_forest$get_species_info())
  clones <- sample_forest$get_species_info() %>% 
    pull(mutant)
  clone_colors <- RColorBrewer::brewer.pal(n_clones, "Dark2")
  names(clone_colors) <- clones
  return(clone_colors)
}

get_karyotypes_colors = function(karyotypes)
{
  karyo_colors = c(
    '1:1' = ggplot2::alpha('seagreen4', .8),
    '1:0' = 'steelblue',
    '0:0' = 'darkblue',
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3',
    '3:0' = 'coral3',
    '3:1' = 'palevioletred',
    '3:2' = 'plum4',
    '3:3' = 'gold2',
    '4:0' = 'lightyellow3',
    '4:1' = 'sandybrown',
    '4:2' = 'tomato2',
    '4:3' = 'darkolivegreen4',
    '4:4' = 'orange3'
  )
  
  missing = setdiff(karyotypes,names(karyo_colors))
  nmissing = length(missing)
  
  
  c(karyo_colors, CNAqc:::nmfy(missing, rep('gray', nmissing)))
}

my_ggplot_theme = function(cex = 1)
{
  cex_opt = getOption('CNAqc_cex', default = 1)
  
  ggplot2::theme_light(base_size = 10 * cex_opt) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(.3 * cex_opt, "cm"),
      panel.background = ggplot2::element_rect(fill = 'white')
    )
}

map_muts_to_karyotype <- function(seq_res,phylo_forest){
  samples <- phylo_forest$get_samples_info()$name
  cna_list <- lapply(samples,function(x){
    c <- phylo_forest$get_bulk_allelic_fragmentation(x) %>% 
      mutate(sample_name=x) %>% 
      mutate(segment_id=paste0(chr,":",begin,":",end))
    
  })
  cna <- do.call("rbind",cna_list)
  muts <- list()
  
  s_seq_long <- seq_res %>% 
    seq_to_long()
  
  for (s in samples){
    mutations <-  s_seq_long %>% 
      filter(sample_name==s)
    cna_sample <- cna %>%
      filter(sample_name==s)
    mapped_mutations <- mutations %>%
      inner_join(cna_sample, by = c("chr","sample_name"), relationship = "many-to-many") %>%  # Allow many-to-many relationships
      filter(from >= begin & to <= end) %>%   # Filter for correct segments
      mutate(segment_id = paste0("chr", chr, ":", begin, ":", end),
             karyotype = paste0(major, ":", minor))
    
    muts[[s]] <-mapped_mutations
  }
  muts_all <- do.call("rbind",muts)
  return(muts_all)
}

absolute_to_relative_coordinates <- function(muts, reference = CNAqc::chr_coordinates_GRCh38){
  vfrom = reference$from
  names(vfrom) = reference$chr
  
  muts %>%
    mutate(
      begin = begin + vfrom[chr],
      end = end + vfrom[chr])
}

blank_genome = function(ref = "GRCh38", 
                        chromosomes = paste0('chr', c(1:22, 'X', 'Y')), 
                        label_chr = -0.5, 
                        cex = 1){
  reference_coordinates = CNAqc::chr_coordinates_GRCh38 %>% filter(chr %in% chromosomes)
  
  low = min(reference_coordinates$from)
  upp = max(reference_coordinates$to)
  
  
  #change the solid and dashed lines for better separating chromosomes.
  p1 = ggplot2::ggplot(reference_coordinates) +
    CNAqc:::my_ggplot_theme(cex = cex) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = centromerStart,
        xend = centromerStart,#centromerEnd,
        y = 0,
        yend = Inf
      ),
      size = .1,
      color = 'black',
      linetype = 8
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = centromerEnd,
        xend = centromerEnd,#centromerEnd,
        y = 0,
        yend = Inf
      ),
      size = .1,
      color = 'black',
      linetype = 8
    )
  
  
  p1 = p1 + ggplot2::geom_rect(
    data = reference_coordinates,
    ggplot2::aes(
      xmin = from,
      xmax = from,
      ymin = 0,
      ymax = Inf
    ),
    alpha = 1,
    colour = 'grey',
  )
  
  p1 = p1 +
    ggplot2::geom_hline(yintercept = 0,
                        size = 1,
                        colour = 'gainsboro') +
    ggplot2::geom_hline(
      yintercept = 1,
      size = .3,
      colour = 'black',
      linetype = 'dashed'
    ) +
    ggplot2::labs(x = "Chromosome",
                  y = "Major/ minor allele") +
    ggpubr::rotate_y_text() +
    # ggpubr::rotate_x_text() +
    # xlim(low, upp) +
    
    #set the chr names in the centromer positions.
    ggplot2::scale_x_continuous(
      breaks = c(0, reference_coordinates$centromerStart, upp),
      labels = c("", gsub(pattern = 'chr', replacement = '', reference_coordinates$chr), "")
    )
  
  return(p1)
}


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


get_exposure_ends <- function(phylo_forest) {
  exposure <- phylo_forest$get_exposures()
  
  time_points <- exposure%>% dplyr::pull(time) %>% unique
  signatures <- exposure %>% dplyr::pull(signature) %>% unique
  for (t in time_points) {
    for (signature in signatures) {
      if (nrow(exposure %>%
               dplyr::filter(.data$time == t,
                             .data$signature == signature)) == 0) {
        exposure[nrow(exposure)+1,] <- c(as.numeric(t), signature,
                                         as.numeric(0), gsub("[0-9]*$","",
                                                             signature))
      }
    }
  }
  
  last_time <- max(phylo_forest$get_samples_info()["time"])
  time_points <- sort(unique(append(time_points, last_time)))
  end_time <- apply(exposure, 1, function(x) {
    time_points[which( round(time_points, 2) == round(as.numeric(x[1]), 2) ) + 1]
  })
  end_time[is.na(end_time)] <- as.numeric(last_time)
  
  exposure$end_time <- end_time
  
  exposure <- exposure[order(exposure$time, exposure$signature), ]
  
  exposure[, c(1, 3)] <- sapply(exposure[, c(1, 3)], as.numeric)
  
  return(exposure)
}

annotate_drivers <- function(phylo_forest){
  ref_path <- phylo_forest$get_reference_path()
  driver_table_path <- gsub(pattern = "reference.fasta",replacement = "drivers.txt",x = ref_path)
  driver_table <-  read.csv(driver_table_path,header=T,sep="\t")
  known_drivers <- driver_table %>%
    dplyr::mutate(chr=as.character(chr)) %>%
    dplyr::rename(driver_label=driver_code) %>% 
    dplyr::rename(pos=from) %>% 
    dplyr::mutate(chr = paste0("chr", chr)) %>% 
    dplyr::select(chr,pos,driver_gene,driver_label)
  
  drivers = phylo_forest$get_driver_mutations() %>%
    dplyr::rename(pos=start) %>%
    dplyr::mutate(chr = paste0("chr", chr)) %>%
    dplyr::left_join(y = known_drivers,by = c("chr","pos"))
  return(drivers)
}

annotate_plots <- function(plot, drivers, ref) {
  reference_genome <- CNAqc:::get_reference(ref = ref)
  vfrom <- reference_genome$from
  names(vfrom) <- reference_genome$chr
  
  if (nrow(drivers) > 0) {
    drivers <- drivers %>%
      dplyr::mutate(pos = pos + vfrom[chr]) %>%
      dplyr::mutate(end = end + vfrom[chr])
    
    driver_SID <- dplyr::filter(drivers, type == "SID")
    driver_CNA <- dplyr::filter(drivers, type == "CNA")
    driver_WGD <- dplyr::filter(drivers, type == "WGD")
    
    # Only call ggplot_build if there's at least one non-empty driver type
    if (nrow(driver_SID) > 0 || nrow(driver_CNA) > 0 || nrow(driver_WGD) > 0) {
      L <- ggplot2::ggplot_build(plot)$layout$panel_params[[1]]
      
      if (nrow(driver_CNA) > 0) {
        driver_CNA$y <- L$y.range[2] * 0.7
        plot <- plot +
          ggplot2::geom_rect(
            data = driver_CNA,
            ggplot2::aes(xmin = pos, xmax = end, fill = CNA_type),
            ymin = 0, ymax = L$y.range[2],
            size = 0.5, alpha = 0.2
          )
      }
      
      if (nrow(driver_SID) > 0) {
        driver_SID$y <- L$y.range[2] * 0.7
        plot <- plot +
          ggplot2::geom_vline(
            data = driver_SID,
            ggplot2::aes(xintercept = pos),
            linetype = 'dashed',
            color = 'black',
            size = 0.3,
            show.legend = FALSE
          ) +
          ggrepel::geom_label_repel(
            data = driver_SID,
            ggplot2::aes(x = pos, y = y, label = driver_label),
            ylim = c(L$y.range[2] * 0.7, NA),
            size = 2,
            nudge_y = 0,
            nudge_x = 0,
            show.legend = FALSE
          )
      }
      
      if (nrow(driver_WGD) > 0) {
        plot <- plot +
          ggplot2::geom_rect(
            data = driver_WGD,
            ggplot2::aes(xmin = L$x.range[1], xmax = L$x.range[2], fill = type),
            ymin = 0, ymax = L$y.range[2],
            size = 0.5, alpha = 0.2
          )
      }
      
      plot <- plot + 
        ggplot2::scale_fill_manual(values = c("A" = "red", "D" = "blue", "WGD" = "grey"))
    }
  }
  
  return(plot)
}


subsample_somatic_mutations <- function(seq_res,fraction){
  # Set desired fraction of total mutations to sample
  fraction_to_sample <- fraction  # 30%
  
  # Compute number of samples per chromosome while keeping proportions
  sample_sizes <- seq_res %>%
    ungroup() %>% 
    count(chr) %>%
    mutate(sample_size = round(n * fraction_to_sample))
  
  # Sample mutations proportionally from each chromosome
  df_sampled <- seq_res %>%
    ungroup() %>%
    group_by(chr) %>%
    group_modify(~ slice_sample(.x, n = sample_sizes$sample_size[sample_sizes$chr == .y$chr])) %>%
    # slice_sample(n = sample_sizes$sample_size[match(chr, sample_sizes$chr)], replace = FALSE) %>%
    ungroup()
  # data <- df_sampled
  df_sampled = df_sampled %>%
    tidyr::pivot_wider(
      names_from = sample_name,
      values_from = c(NV, DP, VAF),
      names_glue = "{sample_name}.{.value}"
    )
  return(df_sampled)
}

label_mutations <- function(model_df){
  mutation_ids <- unique(model_df$mutation_id)
  labels <- c()
  for (mut in mutation_ids){
    muts_data <- model_df %>% 
      # select(mutation_id,SPN01_1.1.VAF,SPN01_1.2.VAF) %>% 
      filter(mutation_id==mut) %>% 
      select(contains("VAF"))
    
    no_zeros <- muts_data[which(muts_data != 0)]
    if (ncol(no_zeros)>1){
      samples <- gsub(pattern = ".VAF",replacement = "",x = colnames(no_zeros))
      label <- paste0(samples, collapse = "_")
      label <- paste0("SHARED_",label)
      # label <- "SHARED"
      labels <- c(labels,label)
    } else if (ncol(no_zeros)==1){
      samples <- gsub(pattern = ".VAF",replacement = "",x = colnames(no_zeros))
      label <- paste0(samples, collapse = "_")
      label <- paste0("PRIVATE_",label)
      labels <- c(labels,label)
    }
  }
  model_df$label <- labels
  return(model_df)
}

get_pairwise_combinations <- function(phylo_forest){
  sample_names <- phylo_forest$get_samples_info()$name
  pairwise_comb <- combn(sample_names, 2, simplify = FALSE)
  return(pairwise_comb)
}

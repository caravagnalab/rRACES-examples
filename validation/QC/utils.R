plot_qc <- function(cnaqc_list, type = 'simple_clonal'){
  if (type == 'simple_clonal'){
    table <- lapply(names(cnaqc_list), function(sample) {
      x <- cnaqc_list[[sample]]
      s <- str_split(sample, '_') %>% unlist()
      
      xqc = CNAqc:::compute_QC_table(x)
      QC_table = xqc$QC_table
      QC_table$sample = paste0(s[2], '_', s[3])
      QC_table$pPASS = xqc$percentage_PASS
      QC_table$NA_tests = xqc$NA_tests
      
      return(QC_table)
    }) %>% bind_rows()
    
    plot = table %>%
      ggplot2::ggplot(aes(x = karyotype, y = sample, fill = paste(QC))) +
      ggplot2::facet_wrap(type ~ .) +
      ggplot2::geom_tile(aes(width = .8, height = .8)) +
      ggplot2::scale_fill_manual(values = c(
        `PASS` = 'seagreen',
        `FAIL` = 'indianred3',
        `NA` = 'gainsboro'
      )) +
      CNAqc:::my_ggplot_theme() +
      ggplot2::labs(x = NULL, y = NULL, title = "Simple clonal CNAs ") +
      ggplot2::guides(fill = ggplot2::guide_legend('QC (NA not available)')) 
    
  } else if (type == 'complex_clonal'){
    table <- lapply(names(cnaqc_list), function(sample) {
      x <- cnaqc_list[[sample]]
      s <- str_split(sample, '_') %>% unlist()
      if (!is.null(x$peaks_analysis$general)){
        QC_table = x$peaks_analysis$general$expected_peak
        QC_table$sample = paste0(s[2], '_', s[3])
        return(QC_table)
      }
    }) %>% bind_rows()
    
    plot <- ggplot() + 
      ggplot2::geom_tile(
        data = table,
        ggplot2::aes(
          y = karyotype,
          x = multiplicity,
          fill = matched,
          width = .8,
          height = .8
        )
      ) +
      CNAqc:::my_ggplot_theme() +
      ggplot2::scale_fill_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'seagreen', `NA` = 'gainsboro')) + 
      ggplot2::labs(x = NULL, y = NULL, title = "Complex clonal CNAs ") +
      facet_grid(.~sample)
  }
  return(plot)
}

get_statistics_qc <- function(cnaqc_list, purity) {
  tt <- lapply(names(cnaqc_list), function(sample) {
    x <- cnaqc_list[[sample]]
    s <- str_split(sample, '_') %>% unlist()
    
    tt = dplyr::tibble(
      sample = paste0(s[2], '_', s[3]),
      n_muts = x$n_mutations, 
      passed_muts = paste0(x$mutations %>% dplyr::filter(QC_PASS == TRUE) %>% nrow(), ' (', (round(x$mutations %>% dplyr::filter(QC_PASS == TRUE) %>% nrow() / x$mutations %>% nrow(), digits = 2))*100, '%)'), 
      #n_failed_muts = x$mutations %>% dplyr::filter(QC_PASS == FALSE) %>% nrow(), 
      n_clonal_cnas = x$n_cna_clonal,
      n_subclonal_cnas = x$n_cna_subclonal, 
      passed_clonal_cnas = paste0(x$cna %>% dplyr::filter(QC_PASS == TRUE) %>% nrow(), ' (', (round(x$cna %>% dplyr::filter(QC_PASS == TRUE) %>% nrow() / x$cna %>% nrow(), digits = 2))*100, '%)'), 
      #failed_clonal_cnas = x$cna %>% dplyr::filter(QC_PASS == FALSE) %>% nrow(), 
      true_purity = purity,
      caller_purity = x$purity, 
      purity_correction = round(x$peaks_analysis$score, 5),
      QC = x$peaks_analysis$QC,
    ) %>% t %>% 
      as.data.frame()
    
    tt = tt %>% 
      dplyr::mutate(info = rownames(tt)) %>% 
      dplyr::rename(values = V1) %>% 
      dplyr::mutate(gg = dplyr::case_when(info == "sample" ~ "sample", 
                                          str_ends(info, "muts") ~ "mutations", 
                                          str_ends(info, "cnas") ~ "cnas", 
                                          str_detect(info, "purity") ~ "purity")) %>% 
      dplyr::mutate(sample = paste0(s[2], '_', s[3]))
    
    return(tt)
  }) %>% dplyr::bind_rows() %>% 
    tidyr::pivot_wider(values_from = values, names_from = sample)
  
  table <- tt %>%
    dplyr::select(gg, info, everything()) %>% 
    dplyr::mutate(info = ifelse(info == "sample", "", info)) %>% 
    dplyr::mutate(gg = ifelse(info == "QC", "QC", gg)) %>% 
    dplyr::mutate(info = ifelse(info == "QC", " ", info)) %>% 
    dplyr::group_by(gg) %>% 
    gt::gt(row_group_as_column = T, process_md = T) %>% 
    gt::tab_header(title = gt::md("**QC summary**")) %>% 
    gt::tab_options(column_labels.hidden = T,  table.layout = "auto")
  
  return(table)
}





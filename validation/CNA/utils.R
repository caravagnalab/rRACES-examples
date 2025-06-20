## Read data
read_ProCESS = function(spn_id,sample_id,coverage,purity){
  
  # Read CNA data
  CNA_races = readRDS(get_process_cna(spn = spn_id,sample = sample_id)) %>%
    mutate(chr = paste0('chr',chr)) %>% dplyr::rename(from=begin,to=end) %>%
    mutate(seg_id = paste0(chr,':',from,':',to)) %>% as_tibble()
  
  # Read SNP data
  mutations = readRDS(get_mutations(spn = spn_id, type = 'tumour', coverage =coverage, purity = purity))
  
  normal = readRDS(get_mutations(spn = spn_id, type = 'normal')) %>% 
    mutate(chr = paste0('chr', chr)) %>% 
    mutate(mut_id = paste(chr, chr_pos, ':')) %>% 
    filter(causes != 'pre-neoplastic')
  
  # Create table of SNPs with associated BAF and DR in absolute coordinates
  snps = mutations %>% 
    filter(classes == "germinal", causes=="") %>% 
    mutate(chr = paste0('chr', chr)) %>% 
    mutate(mut_id = paste(chr, chr_pos, ':'))
  
  snps = snps %>% ungroup() %>% sample_n(size = nrow(snps)/100)
  normal = normal %>% ungroup() %>% filter(mut_id %in% snps$mut_id)
  
  snps = snps %>%
    rename(pos=chr_pos) %>% select(c("chr", "pos", 
                                     paste0(sample_id,".occurrences"), 
                                     paste0(sample_id,".coverage"), 
                                     paste0(sample_id,".VAF")),
                                   "mut_id")
  colnames(snps)=c("chr", "pos", "NV", "DP", "VAF", "mut_id")
  
  normal = normal  %>% 
    rename(pos=chr_pos) %>% 
    select(c("chr", "pos","normal_sample.occurrences","normal_sample.coverage", "normal_sample.VAF", "mut_id"))
  colnames(normal)=c("chr", "pos", "NV_normal", "DP_normal", "VAF_normal", "mut_id")
  
  
  joint_table_snps = left_join(snps, normal, by=c('chr', 'pos'))
  coverage_normal = 30
  coverage_tumor = as.integer(coverage)
  norm_const = coverage_normal/coverage_tumor
  joint_table_snps = joint_table_snps %>% 
    rename(BAF = VAF) %>% 
    mutate(DR = (DP/DP_normal)*norm_const)
  joint_table_snps_shifted=lapply(chromosomes, function(c){
    from_c = CNAqc:::get_reference('hg38') %>% filter(chr == c) %>% pull(from)
    joint_table_snps %>% filter(chr == c) %>% mutate(pos = pos+from_c)
  })
  joint_table_snps_shifted = Reduce(rbind,joint_table_snps_shifted)
  
  return(list('CNA'=CNA_races, 'snps'=joint_table_snps_shifted))
}

read_ASCAT = function(spn_id,sample_id,coverage,purity){
  caller = 'ascat'
  ascat_results <- get_sarek_cna_file(spn = spn_id,
                                      coverage = coverage,
                                      purity = purity,
                                      caller = caller,
                                      type = "tumour",
                                      sampleID = sample_id)
  
  CNA_ascat = read.csv(ascat_results$cnvs, sep='\t') %>%
    mutate(chr = paste0('chr',chr)) %>% rename(from=startpos,to=endpos,major=nMajor,minor=nMinor) %>% as_tibble()
  # BAF and DR ascat
  BAF_file = data.table::fread(ascat_results$tumourBAF, sep='\t')
  colnames(BAF_file) = c('id', 'Chromosome', 'Position', 'BAF')
  DR_file = data.table::fread(ascat_results$tumourLogR, sep='\t')
  colnames(DR_file) = c('id', 'Chromosome', 'Position', 'LogR')
  purity_ploidy = read.csv(ascat_results$purityploidy, sep='\t')
  
  # Join and shif BAF and DR
  BAF_DR = left_join(
    BAF_file,
    DR_file %>% mutate(DR = exp(LogR))) %>%
    mutate(Chromosome = paste0('chr', Chromosome))
  random_indeces = sample(BAF_DR$Position, as.integer(length(BAF_DR$Position)/100), replace = FALSE)
  BAF_DR_filtered = BAF_DR %>% filter(Position %in% random_indeces)
  BAF_DR_shifted = lapply(1:nrow(BAF_DR_filtered), function(r){
    #print(r)
    chromosome = BAF_DR_filtered[r,]$Chromosome
    from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
    BAF_DR_filtered[r,] %>% mutate(Position = Position + from_chromosome)
  })
  BAF_DR_shifted = Reduce(rbind, BAF_DR_shifted)
  
  return(list('CNA'=CNA_ascat, 'BAF_DR'=BAF_DR_shifted %>% as_tibble(), 'purity_ploidy'= purity_ploidy))
}

read_Sequenza = function(spn_id,sample_id,coverage,purity, seq_dir){
  caller = 'sequenza'
  sequenza_results <- get_sarek_cna_file(spn = spn_id,
                                         coverage = coverage,
                                         purity = purity,
                                         caller = caller,
                                         type = "tumour",
                                         sampleID = sample_id)
  CNA_sequenza = read.csv(sequenza_results$segments, sep='\t')
  CNA_sequenza = CNA_sequenza %>% select(!c(N.BAF,sd.BAF,N.ratio,sd.ratio,CNt,LPP))
  colnames(CNA_sequenza) = c('chr', 'from', 'to','BAF','DR','major','minor')
  BAF_DR_file_sequenza = sequenza::read.seqz(paste0(seq_dir,'/filtered.seqz'))
  BAF_DR_file_sequenza_shifted = lapply(chromosomes, function(c){
    from_c = CNAqc:::get_reference('hg38') %>% filter(chr == c) %>% pull(from)
    BAF_DR_file_sequenza %>% filter(chromosome == c) %>% 
      mutate(shifted_position = position+from_c)
  })
  BAF_DR_file_sequenza_shifted = Reduce(rbind, BAF_DR_file_sequenza_shifted)
  purity_ploidy_sequenza = read.csv(sequenza_results$confints_CP, sep='\t')
  
  return(list('CNA'=CNA_sequenza, 'BAF_DR'=BAF_DR_file_sequenza_shifted, 'purity_ploidy'= purity_ploidy_sequenza))
}

read_CNVkit = function(spn_id,sample_id,coverage,purity){
  caller='cnvkit'
  cnvkit_results <- get_sarek_cna_file(spn = spn_id,
                                       coverage = coverage,
                                       purity = purity,
                                       caller = caller,
                                       type = "tumour",
                                       sampleID = sample_id)
  
  CNA_cnvkit = read.csv(cnvkit_results$somatic.call,sep='\t')
  CNA_cnvkit = CNA_cnvkit %>% select(!c(gene,ci_hi,ci_lo,depth,probes,weight))
  colnames(CNA_cnvkit) = c('chr', 'from','to','logR','CN')
  # DR cnvkit
  DR_file_cnvkit = data.table::fread(cnvkit_results$cnr, sep='\t') %>% as_tibble()
  
  # ## inspect CNVkit calls
  # CNA_cnvkit_shifted =lapply(chromosomes, function(c){
  #   from_c = CNAqc:::get_reference('hg38') %>% filter(chr == c) %>% pull(from)
  #   CNA_cnvkit %>% filter(chromosome == c) %>% 
  #     mutate(start = start+from_c,
  #            end = end+from_c)
  # })
  # CNA_cnvkit_shifted = Reduce(rbind,CNA_cnvkit_shifted)
  
  DR_file_cnvkit_shifted =lapply(chromosomes, function(c){
    from_c = CNAqc:::get_reference('hg38') %>% filter(chr == c) %>% pull(from)
    DR_file_cnvkit %>% filter(chromosome == c) %>% 
      mutate(start = start+from_c,
             end = end+from_c)
  })
  DR_file_cnvkit_shifted = Reduce(rbind,DR_file_cnvkit_shifted)
  
  return(list('CNA'=CNA_cnvkit, 'DR'=DR_file_cnvkit_shifted))
}

## Create joint segmentation
create_joint_segmentation = function(CNA_ProCESS, CNA_target, caller, chromosomes){
  
  joint_segmentation = lapply(chromosomes, function(c){
    #print(c)
    CNA_ProCESS_chr = CNA_ProCESS %>% filter(chr==c) %>% 
      mutate(group = cumsum(
        lag(major, default = dplyr::first(major)) != major |
          lag(minor, default = dplyr::first(minor)) != minor |
          lag(ratio, default = dplyr::first(ratio)) != ratio
      )) %>%
      group_by(chr, major, minor, ratio, group) %>%
      summarise(from = min(from), to = max(to), .groups = "drop") %>%
      select(chr, from, to, major, minor, ratio) %>% 
      arrange(from, to)
    
    CNA_target_chr = CNA_target %>% filter(chr==c)
    froms = c(CNA_ProCESS_chr$from %>% unique(), CNA_target_chr$from %>% unique())
    tos = c(CNA_ProCESS_chr$to %>% unique(), CNA_target_chr$to %>% unique())
    breakpoints = c(froms, tos) %>% sort()
    df = data.frame()
    
    for (i in 1:(length(breakpoints)-1)){
      f = breakpoints[i]
      t = breakpoints[i+1]
      true_cna= CNA_ProCESS_chr %>% filter(from <= f, from < t, to>=t) 
      inferred_cna= CNA_target_chr %>% filter(from <= f, from < t, to>=t) 
      if (nrow(true_cna)==0){
        true_M1=NA
        true_m1=NA
        true_M2=NA
        true_m2=NA
        true_cn=NA
      }
      if (nrow(true_cna)==1){
        true_M1=true_cna %>% pull(major)
        true_m1=true_cna %>% pull(minor)
        true_M2=NA
        true_m2=NA
        true_cn=true_M1+true_m1
      }
      if (nrow(true_cna)>1){
        true_M1=true_cna[1,] %>% pull(major)
        true_m1=true_cna[1,] %>% pull(minor)
        true_M2=true_cna[2,] %>% pull(major)
        true_m2=true_cna[2,] %>% pull(minor)
        ratio=true_cna[1,] %>% pull(ratio)
        true_cn = ratio*(true_M1+true_m1) + (1-ratio)*(true_M2+true_m2) 
      }
      if (nrow(inferred_cna)==0){
        inferred_M=NA
        inferred_m=NA
        inferred_cn=NA
      }
      if (nrow(inferred_cna)==1){
        if (caller!='cnvkit'){
          inferred_M=inferred_cna %>% pull(major)
          inferred_m=inferred_cna %>% pull(minor)
          inferred_cn=inferred_M+inferred_m
        }else{
          inferred_M=NA
          inferred_m=NA
          inferred_cn=inferred_cna %>% pull(CN)
        }
      }
      df = rbind(df, data.frame(
        'from' =f,
        'to'=t,
        'chromosome'=c,
        'TRUE_Major1'= true_M1,
        'TRUE_minor1'= true_m1,
        'TRUE_Major2'= true_M2,
        'TRUE_minor2'= true_m2,
        'TRUE_CN' = true_cn,
        'INFERRED_Major1' = inferred_M,
        'INFERRED_minor1' = inferred_m,
        'INFERRED_CN' = inferred_cn,
        'caller'=caller
      ))
    }
    df
    
  }) #%>% Reduce(rbind)
  joint_segmentation = Reduce(rbind, joint_segmentation)
  
  joint_segmentation_shifted = lapply(1:nrow(joint_segmentation), function(r){
    chromosome = joint_segmentation[r,]$chromosome
    from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
    joint_segmentation[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
  })
  joint_segmentation_shifted = Reduce(rbind, joint_segmentation_shifted)
  
  if (caller %in% c('ascat','sequenza')){
    joint_segmentation_shifted_longer = joint_segmentation_shifted %>% 
      mutate(is_match = 
               case_when(
                 TRUE_Major1==INFERRED_Major1 & TRUE_minor1==INFERRED_minor1 & is.na(TRUE_Major2) | 
                   TRUE_Major1==INFERRED_Major1 & TRUE_minor1==INFERRED_minor1 & TRUE_Major2==INFERRED_Major1 & TRUE_minor2==INFERRED_minor1 ~ 'complete match',
                 TRUE_Major1==INFERRED_Major1 & TRUE_minor1==INFERRED_minor1 & !is.na(TRUE_Major2) ~ 'undetected subclone',
                 .default = 'no match'
               )) %>%
      pivot_longer(
        cols = c('TRUE_Major1', 'TRUE_minor1', 'TRUE_Major2', 'TRUE_minor2', 'INFERRED_Major1', 'INFERRED_minor1'),  
        names_to = "Type",  # New column for variable names
        values_to = "Value"  # New column for values
      )
  }else{
    joint_segmentation_shifted_longer = joint_segmentation_shifted %>% 
      mutate(is_match = 
               case_when(
                 abs(TRUE_CN - INFERRED_CN)<.1 ~ 'complete match',
                 abs(TRUE_CN - INFERRED_CN)>=.1 & abs(TRUE_CN - INFERRED_CN)<.5 ~ 'close',
                 .default = 'no match'
               )) %>% select(from, to, chromosome, TRUE_CN, INFERRED_CN,is_match)%>%
      pivot_longer(
        cols = c('TRUE_CN', 'INFERRED_CN'),  
        names_to = "Type",  # New column for variable names
        values_to = "Value"  # New column for values
      )
  }
  
  return(list('joint_segmentation'=joint_segmentation_shifted,
              'joint_segmentation_long'=joint_segmentation_shifted_longer))
}

## Shift segments
shift_segments = function(CNA){
  CNA_shifted = lapply(1:nrow(CNA), function(r){
  chromosome = CNA[r,]$chr
  from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
  
  if ("original_from" %in% colnames(CNA)){
    CNA[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome,
                       original_from = original_from + from_chromosome, original_to = original_to + from_chromosome)
  }else{
    CNA[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
  }
  
  })
  CNA_shifted = Reduce(rbind, CNA_shifted)
  CNA_shifted
}

## Compute correctness 
compute_correctness = function(df){
  1-( (df %>% filter(is_match == 'no match') %>% mutate(len=to-from) %>%
         pull(len) %>% unique() %>% sum()) / (df %>% mutate(len=to-from) %>%
                                                pull(len) %>% unique() %>% sum()))
}

## Expected BAF and DR functions
expected_BAF = function(nA1, nB1, nA2, nB2, purity, ccf){
  if (is.na(nA2) & is.na(nB2)){
    nA2 = 0
    nB2 = 0
  }
  nom = min(nA1*ccf + nA2*(1-ccf), nB1*ccf + nB2*(1-ccf))*purity + (1-purity)
  den = ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*purity + 2*(1-purity)
  nom/den
}
expected_DR = function(nA1, nB1, nA2, nB2, purity, ccf, ploidy){
  if (is.na(nA2) & is.na(nB2)){
    nA2 = 0
    nB2 = 0
  }
  nom = ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*purity + 2*(1-purity)
  nom/(ploidy*purity + 2*(1-purity))
}

## True Ploidy
compute_true_ploidy = function(CNA_ProCESS){
  # Compute real ploidy
  CNA_ProCESS = CNA_ProCESS %>% mutate(len=to-from)
  genome_len =CNA_ProCESS$len %>% unique() %>% sum()
  CNA_ProCESS = CNA_ProCESS %>% mutate(seg_id = paste0(chr, ':', from,':',to))
  segments = CNA_ProCESS$seg_id %>% unique()
  ploidy = 0
  # purity_number = as.double(purity)
  for (s in segments){
    seg = CNA_ProCESS %>% filter(seg_id == s)
    nA1 = seg[1,]$major
    nB1 = seg[1,]$minor
    genome_fraction = seg[1,]$len/genome_len
    if (nrow(seg)>1){
      nA2 = seg[2,]$major
      nB2 = seg[2,]$minor
      ccf = seg[1,]$ratio
    }else{
      nA2 = 0
      nB2 = 0
      ccf = 1
    }
    
    # ploidy = ploidy + ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*purity_number*genome_fraction
    ploidy = ploidy + ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*genome_fraction
  }
  # ploidy = ploidy + 2*(1-purity_number)
  ploidy
}

## Theoretical vs inferred BAF and DR (ascat)
theoretical_vs_inferred_BAF_DR = function(CNA_ProCESS,CNA_ascat,purity_number,ploidy,purity_ploidy){
  already_seen_segment = c()
  CNA_ProCESS = lapply(1:nrow(CNA_ProCESS), function(r){
    chromosome = CNA_ProCESS[r,]$chr
    from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
    CNA_ProCESS[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
  })
  CNA_ProCESS = Reduce(rbind, CNA_ProCESS)
  transposed_CNA_ProCESS= lapply(CNA_ProCESS$seg_id%>%unique(), function(s){
    CNA_ProCESS_s = CNA_ProCESS %>% filter(seg_id==s)
    if (nrow(CNA_ProCESS_s)>1){
      data.frame("chr"=CNA_ProCESS_s[1,]$chr,
                 "from"=CNA_ProCESS_s[1,]$from,
                 "to"=CNA_ProCESS_s[1,]$to,
                 "Major1"=CNA_ProCESS_s[1,]$major,
                 "minor1"=CNA_ProCESS_s[1,]$minor,
                 "Major2"=CNA_ProCESS_s[2,]$major,
                 "minor2"=CNA_ProCESS_s[2,]$minor,
                 "ratio"=CNA_ProCESS_s[1,]$ratio)
    }else{
      data.frame("chr"=CNA_ProCESS_s[1,]$chr,
                 "from"=CNA_ProCESS_s[1,]$from,
                 "to"=CNA_ProCESS_s[1,]$to,
                 "Major1"=CNA_ProCESS_s[1,]$major,
                 "minor1"=CNA_ProCESS_s[1,]$minor,
                 "Major2"=NA,
                 "minor2"=NA,
                 "ratio"=CNA_ProCESS_s[1,]$ratio)
    }
  })
  transposed_CNA_ProCESS = Reduce(rbind,transposed_CNA_ProCESS)
  CNA_ascat = lapply(1:nrow(CNA_ascat), function(r){
    chromosome = CNA_ascat[r,]$chr
    from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
    CNA_ascat[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
  })
  CNA_ascat = Reduce(rbind, CNA_ascat)
  
  simulated_baf_dr = lapply(1:nrow(transposed_CNA_ProCESS), function(r){
    simulated_baf=expected_BAF(transposed_CNA_ProCESS[r,]$Major1, transposed_CNA_ProCESS[r,]$minor1, transposed_CNA_ProCESS[r,]$Major2, transposed_CNA_ProCESS[r,]$minor2, purity_number, transposed_CNA_ProCESS[r,]$ratio)
    simulated_dr=expected_DR(transposed_CNA_ProCESS[r,]$Major1, transposed_CNA_ProCESS[r,]$minor1,transposed_CNA_ProCESS[r,]$Major2, transposed_CNA_ProCESS[r,]$minor2, purity_number, transposed_CNA_ProCESS[r,]$ratio, ploidy)
    cbind(transposed_CNA_ProCESS[r,], data.frame('simulated_baf'=c(simulated_baf), 'simulated_dr'=c(simulated_dr)))
  })
  simulated_baf_dr = Reduce(rbind, simulated_baf_dr)
  
  inferred_baf_dr = lapply(1:nrow(CNA_ascat), function(r){
    inferred_baf=expected_BAF(CNA_ascat[r,]$major, CNA_ascat[r,]$minor, 0, 0, purity_ploidy$AberrantCellFraction, 1)
    inferred_dr=expected_DR(CNA_ascat[r,]$major, CNA_ascat[r,]$minor, 0, 0, purity_ploidy$AberrantCellFraction, 1, purity_ploidy$Ploidy)
    cbind(CNA_ascat[r,], data.frame('inferred_baf'=c(inferred_baf), 'inferred_dr'=c(inferred_dr)))
  })
  inferred_baf_dr = Reduce(rbind,inferred_baf_dr)
  
  return(list('simulated_baf_dr' = simulated_baf_dr,'inferred_baf_dr_ascat' = inferred_baf_dr))
}

## Compute segmentation correctness
covered_genome = function(CNA_target, chromosome){
  chr_len = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(length)
  covered = CNA_target %>% filter(chr == chromosome) %>% mutate(l=to-from) %>% pull(l) %>% sum()
  (covered / chr_len)*100
}
breakpoint_analysis = function(CNA_ProCESS, CNA_target, chromosome, th=1e7){
  ProCESS_BP = CNA_ProCESS %>% filter(chr==chromosome) #%>% select(from, to)
  BP_f = ProCESS_BP$from
  BP_t = ProCESS_BP$to
  
  CNA_target_c = CNA_target %>% filter(chr==chromosome)
  BP_target_f = CNA_target_c$from %>% unique()
  BP_target_t = CNA_target_c$to %>% unique()
  
  missed_bp = 0
  missed_bp_coord = c()
  added_bp = 0
  added_bp_coord = c()
  distances = c()
  matching_coord = c()
  
  for (bp in BP_f){
    dist = min(abs(bp - BP_target_f))
    if (dist > th){
      missed_bp = missed_bp+1
      missed_bp_coord = c(missed_bp_coord,bp)
    }else{
        distances=c(distances,dist)
        matching_coord = c(matching_coord, bp)
        }
  }
  
  for (bp in BP_t){
    dist = min(abs(bp - BP_target_t))
    if (dist > th){
      missed_bp = missed_bp+1
      missed_bp_coord = c(missed_bp_coord,bp)
    }else{
        distances=c(distances,dist)
        matching_coord = c(matching_coord, bp)
        }
  }
  
  for (bp in BP_target_f){
    dist = min(abs(bp - BP_f))
    if (dist > th){
      added_bp = added_bp+1
      added_bp_coord = c(added_bp_coord, bp)
      }
  }
  
  list(
    'chromosome' = chromosome,
    'missed_bp_coord' = missed_bp_coord %>% unique(),
    'added_bp_coord' = added_bp_coord %>% unique(),
    'n_missed' = missed_bp,
    'n_added' = added_bp,
    'distances_between_matching' = distances %>% unique(),
    'matching_coord' = matching_coord %>% unique()
  )
}
segmentation_analysis = function(CNA_ProCESS, CNA_target, chromosomes, th=1e7){
  bp_df = lapply(chromosomes, function(c){
    breakpoint_analysis(CNA_ProCESS, CNA_target, c, th)
  })
}


## Find matching segments between two different segmentations
# find_matching_segments = function(CNA_ProCESS_unique,CNA_target,seg_ids){
#   best_segments = lapply(seg_ids, function(s){
#     # print(s)
#     f = CNA_ProCESS_unique %>% filter(seg_id==s) %>% pull(from)
#     t = CNA_ProCESS_unique %>% filter(seg_id==s) %>% pull(to)
#     c = CNA_ProCESS_unique %>% filter(seg_id==s) %>% pull(chr)
#     max_overlap_df = CNA_target %>% filter(chr == c) %>% rowwise() %>% mutate(l= (min(t,to)-max(f,from))) %>%
#       #mutate(seg_len=to-from) %>% rowwise() %>%
#       mutate(p= min( (l/(t-f))*100, 100) )
#     max_overlap = max(max_overlap_df$l)
#     best_match = max_overlap_df %>% filter(l==max_overlap) 
#     if (nrow(best_match)>0){
#       best_match=best_match %>% mutate(original_from= f, original_to= t, distance = mean(abs(from-f), abs(to-t)))
#     }
#   })
#   best_segments = Reduce(rbind,best_segments)
#   best_segments = best_segments %>% rowwise() %>% mutate(ids = paste0(chr, ':',from,':',to))
#   ids = best_segments$ids %>% unique()
#   filtered_best_segments = lapply(ids, function(i){
#     best_segments_i= best_segments %>% filter(ids == i)
#     best_segments_i = best_segments_i %>% filter(p==max(best_segments_i$p)) 
#     # best_segments_i = best_segments_i %>% filter(l==max(best_segments_i$l))
#     best_segments_i
#   })
#   filtered_best_segments = Reduce(rbind, filtered_best_segments)
#   filtered_best_segments
# }







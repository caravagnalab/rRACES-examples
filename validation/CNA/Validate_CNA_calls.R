options(bitmapType='cairo')
library(dplyr)
library(ProCESS)
library(optparse)
library(tidyr)
library(ggplot2)

############ Parse command-line arguments
option_list <- list(make_option(c("--sample_id"), type = "character", default = 'SPN03_1.1'),
		                make_option(c("--spn_id"), type = "character", default = 'SPN03'),
                    make_option(c("--purity"), type = "character", default = '0.9'),
                    make_option(c("--coverage"), type = "character", default = '100'),
                    make_option(c("--purity_th"), type = "character", default = '.1'),
                    make_option(c("--correct_th"), type = "character", default = '.6'))
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'

sample_id = opt$sample_id
spn_id = opt$spn_id
coverage = paste0(opt$coverage,"x")
purity = paste0(opt$purity,"p")
purity_th = opt$purity_th
correct_th = opt$correct_th
caller = "ascat"


gender_file <- read.table(file = paste0(data_dir,spn_id,"/process/subject_gender.txt"),header = FALSE,col.names = "gender")
gender <- gender_file$gender

if (gender=="XX"){
  chromosomes = c(paste0('chr',1:22), 'chrX')
} else {
  chromosomes = c(paste0('chr',1:22), 'chrX', 'chrY')
}

############ Load data
#### ProCESS data
# CNA calls
CNA_races = readRDS(paste0(data_dir,spn_id,'/process/cna_data/',sample_id,'_cna.rds')) %>%
  mutate(chr = paste0('chr',chr)) %>% dplyr::rename(from=begin,to=end) %>% as_tibble()

# SNP data from which compute BAF and DR
mutations = readRDS(paste0(data_dir, '/',spn_id,'/sequencing/tumour/purity_',strsplit(purity,'p')[[1]],
                           '/data/mutations/seq_results_muts_merged_coverage_',coverage,'.rds'))

normal = readRDS(paste0(data_dir, '/',spn_id,'/sequencing/normal/purity_1', 
                        '/data/mutations/seq_results_muts_merged_coverage_30x.rds')) %>% 
  mutate(chr = paste0('chr', chr)) %>% 
  mutate(mut_id = paste(chr, chr_pos, ':')) %>% 
  filter(causes != 'pre-neoplastic')

message("Reading ProCESS data")
#### ascat data
# CNA calls
CNA_ascat = read.csv(paste0(data_dir, spn_id, '/sarek/', coverage, '_', purity, 
                      '/variant_calling/ascat/',sample_id,'_vs_normal_sample/',
                      sample_id,'_vs_normal_sample.cnvs.txt'), sep='\t') %>%
  mutate(chr = paste0('chr',chr)) %>% rename(from=startpos,to=endpos,major=nMajor,minor=nMinor) %>% as_tibble()
# BAF and DR ascat
BAF_file = data.table::fread(paste0(data_dir, spn_id, '/sarek/', coverage, '_', purity, 
                      '/variant_calling/ascat/',sample_id,'_vs_normal_sample/',
                      sample_id,'_vs_normal_sample.tumour_tumourBAF.txt'), sep='\t')
colnames(BAF_file) = c('id', 'Chromosome', 'Position', 'BAF')
DR_file = data.table::fread(paste0(data_dir, spn_id, '/sarek/', coverage, '_', purity, 
                  '/variant_calling/ascat/',sample_id,'_vs_normal_sample/',
                  sample_id,'_vs_normal_sample.tumour_tumourLogR.txt'), sep='\t')
colnames(DR_file) = c('id', 'Chromosome', 'Position', 'LogR')
purity_ploidy = read.csv(paste0(data_dir, spn_id, '/sarek/', coverage, '_', purity, 
                                '/variant_calling/ascat/',sample_id,'_vs_normal_sample/',
                                sample_id,'_vs_normal_sample.purityploidy.txt'), sep='\t')
message("Reading ASCAT data")
#### CNVkit data

############ Process data

message("Create joint ProCESS SNPs table")
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
coverage_tumor = as.integer(strsplit(coverage, 'x')[[1]])
norm_const = coverage_normal/coverage_tumor
joint_table_snps = joint_table_snps %>% 
  rename(BAF = VAF) %>% 
  mutate(DR = (DP/DP_normal)*norm_const)
joint_table_snps_shifted=lapply(chromosomes, function(c){
  from_c = CNAqc:::get_reference('hg38') %>% filter(chr == c) %>% pull(from)
  joint_table_snps %>% filter(chr == c) %>% mutate(pos = pos+from_c)
})
joint_table_snps_shifted = Reduce(rbind,joint_table_snps_shifted)


message("Create joint table ProCESS and ascat calls") 
joint_segmentation = lapply(chromosomes, function(c){
  print(c)
  CNA_races_chr = CNA_races %>% filter(chr==c) %>% 
    mutate(group = cumsum(
      lag(major, default = dplyr::first(major)) != major |
        lag(minor, default = dplyr::first(minor)) != minor |
        lag(ratio, default = dplyr::first(ratio)) != ratio
    )) %>%
    group_by(chr, major, minor, ratio, group) %>%
    summarise(from = min(from), to = max(to), .groups = "drop") %>%
    select(chr, from, to, major, minor, ratio) %>% 
    arrange(from, to)

  CNA_ascat_chr = CNA_ascat %>% filter(chr==c)
  froms = c(CNA_races_chr$from %>% unique(), CNA_ascat_chr$from %>% unique())
  tos = c(CNA_races_chr$to %>% unique(), CNA_ascat_chr$to %>% unique())
  breakpoints = c(froms, tos) %>% sort()
  df = data.frame()
  for (i in 1:(length(breakpoints)-1)){
    f = breakpoints[i]
    t = breakpoints[i+1]
    true_cna= CNA_races_chr %>% filter(from <= f, from < t, to>=t) 
    inferred_cna= CNA_ascat_chr %>% filter(from <= f, from < t, to>=t) 
    if (nrow(true_cna)==0){
      true_M1=NA
      true_m1=NA
      true_M2=NA
      true_m2=NA
    }
    if (nrow(true_cna)==1){
      true_M1=true_cna %>% pull(major)
      true_m1=true_cna %>% pull(minor)
      true_M2=NA
      true_m2=NA
    }
    if (nrow(true_cna)>1){
      true_M1=true_cna[1,] %>% pull(major)
      true_m1=true_cna[1,] %>% pull(major)
      true_M2=true_cna[2,] %>% pull(major)
      true_m2=true_cna[2,] %>% pull(major)
    }
    if (nrow(inferred_cna)==0){
      true_M=NA
      true_m=NA
    }
    if (nrow(inferred_cna)==1){
      true_M=inferred_cna %>% pull(major)
      true_m=inferred_cna %>% pull(minor)
    }
    df = rbind(df, data.frame(
      'from' =f,
      'to'=t,
      'chromosome'=c,
      'TRUE_Major1'= true_M1,
      'TRUE_minor1'= true_m1,
      'TRUE_Major2'= true_M2,
      'TRUE_minor2'= true_m2,
      'INFERRED_Major1'=true_M,
      'INFERRED_minor1'=true_m
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

# Compute CNA correctness - percentage of the genome correctly called by ascat 
# (call is considered correct if the major and minor allele match, event if the match 
# is reffered only to the most prevalent subclone)
message("Compute metrics")
ascat_correctness =  1-( (joint_segmentation_shifted_longer %>% filter(is_match == 'no match') %>% mutate(len=to-from) %>%
  pull(len) %>% unique() %>% sum()) / (joint_segmentation_shifted_longer %>% mutate(len=to-from) %>%
                                         pull(len) %>% unique() %>% sum()))
purity_correctness = purity_ploidy$AberrantCellFraction - as.double(strsplit(purity, 'p')[[1]])
if (purity_th < purity_correctness & correct_th < ascat_correctness){state = 'PASS'}else{state='FAIL'}

## BAF and DR ascat
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
  nom/ploidy
}

# Compute real ploidy
CNA_races = CNA_races %>% mutate(len=to-from)
genome_len =CNA_races$len %>% unique() %>% sum()
CNA_races = CNA_races %>% mutate(seg_id = paste0(chr, ':', from,':',to))
segments = CNA_races$seg_id %>% unique()
ploidy = 0
purity_number = as.double(strsplit(purity, 'p')[[1]])
for (s in segments){
  seg = CNA_races %>% filter(seg_id == s)
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
  
  ploidy = ploidy + ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*purity_number*genome_fraction
}
# ploidy = ploidy + 2*(1-purity_number)

# Theoretical vs inferred BAF and DR
already_seen_segment = c()
CNA_races = lapply(1:nrow(CNA_races), function(r){
  chromosome = CNA_races[r,]$chr
  from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
  CNA_races[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
})
CNA_races = Reduce(rbind, CNA_races)
transposed_CNAraces= lapply(CNA_races$seg_id%>%unique(), function(s){
  CNA_races_s = CNA_races %>% filter(seg_id==s)
  if (nrow(CNA_races_s)>1){
    data.frame("chr"=CNA_races_s[1,]$chr,
               "from"=CNA_races_s[1,]$from,
               "to"=CNA_races_s[1,]$to,
               "Major1"=CNA_races_s[1,]$major,
               "minor1"=CNA_races_s[1,]$minor,
               "Major2"=CNA_races_s[2,]$major,
               "minor2"=CNA_races_s[2,]$minor,
               "ratio"=CNA_races_s[1,]$ratio)
  }else{
    data.frame("chr"=CNA_races_s[1,]$chr,
               "from"=CNA_races_s[1,]$from,
               "to"=CNA_races_s[1,]$to,
               "Major1"=CNA_races_s[1,]$major,
               "minor1"=CNA_races_s[1,]$minor,
               "Major2"=NA,
               "minor2"=NA,
               "ratio"=CNA_races_s[1,]$ratio)
  }
})
transposed_CNAraces = Reduce(rbind,transposed_CNAraces)
CNA_ascat = lapply(1:nrow(CNA_ascat), function(r){
  chromosome = CNA_ascat[r,]$chr
  from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
  CNA_ascat[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
})
CNA_ascat = Reduce(rbind, CNA_ascat)

simulated_baf_dr = lapply(1:nrow(transposed_CNAraces), function(r){
  simulated_baf=expected_BAF(transposed_CNAraces[r,]$Major1, transposed_CNAraces[r,]$minor1, transposed_CNAraces[r,]$Major2, transposed_CNAraces[r,]$minor2, purity_number, transposed_CNAraces[r,]$ratio)
  simulated_dr=expected_DR(transposed_CNAraces[r,]$Major1, transposed_CNAraces[r,]$minor1,transposed_CNAraces[r,]$Major2, transposed_CNAraces[r,]$minor2, purity_number, transposed_CNAraces[r,]$ratio, ploidy)
  cbind(transposed_CNAraces[r,], data.frame('simulated_baf'=c(simulated_baf), 'simulated_dr'=c(simulated_dr)))
})

simulated_baf_dr = Reduce(rbind, simulated_baf_dr)

inferred_baf_dr = lapply(1:nrow(CNA_ascat), function(r){
  inferred_baf=expected_BAF(CNA_ascat[r,]$major, CNA_ascat[r,]$minor, 0, 0, purity_ploidy$AberrantCellFraction, 1)
  inferred_dr=expected_DR(CNA_ascat[r,]$major, CNA_ascat[r,]$minor, 0, 0, purity_ploidy$AberrantCellFraction, 1, purity_ploidy$Ploidy)
  cbind(CNA_ascat[r,], data.frame('inferred_baf'=c(inferred_baf), 'inferred_dr'=c(inferred_dr)))
})
inferred_baf_dr = Reduce(rbind,inferred_baf_dr)

########### Plots
message("Generate plots")
color_by_state = c("TRUE_Major1"=alpha('firebrick', 1),
                   "TRUE_minor1"=alpha('#000080ff',1), 
                   "TRUE_Major2"=alpha('#ff00abb3'),
                   "TRUE_minor2"=alpha('slateblue', 1),
                   "INFERRED_Major1"=alpha('firebrick', .5),
                   "INFERRED_minor1"=alpha('#000080ff',.5))

fill_by_match = c('complete match'= alpha('gainsboro', .03),
                  'undetected subclone'= alpha('goldenrod', .08),
                  'no match' = alpha('indianred', .08)
)
# fill_by_match_2 = c('complete match'= alpha('forestgreen', .5),
#                     'undetected subclone'= alpha('#bcff5cb3', .5),
#                     'no match' = alpha('indianred', .5)
# )
# guides_colors = fill_by_match_2[joint_segmentation_shifted_longer$is_match %>% unique()]


#CNAqc:::blank_genome() + geom_point(data=snps, aes(x=chr_pos, y =SPN01_1.1.VAF))

baf_ascat = CNAqc:::blank_genome() + 
  geom_point(data = BAF_DR_shifted %>% sample_n(size = n()/2), #%>% mutate(BAF = ifelse(BAF>.5, 1-BAF,BAF))
             aes(x=Position, y=BAF), size = .1, alpha=1)+
  geom_segment(data = inferred_baf_dr, aes(x=from, xend=to, y=inferred_baf), color='red', size=1)+
  #geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf), color='steelblue', size=1)+
  #ylim(0,.5)+
  ylab('BAF ascat')

dr_ascat = CNAqc:::blank_genome() + 
  geom_point(data = BAF_DR_shifted %>% sample_n(size = n()/2), 
             aes(x=Position, y=DR), size = .1, alpha=1)+
  geom_segment(data = inferred_baf_dr, aes(x=from, xend=to, y=inferred_dr), color='red', size=1)+
  ylab('DR ascat')+
  ylim(0,4)

baf_races = CNAqc:::blank_genome() + 
  geom_point(data = joint_table_snps_shifted %>% sample_n(size = n()/5), #%>% mutate(BAF = ifelse(BAF>.5, 1-BAF,BAF)) 
             aes(x=pos, y=BAF), size = .1, alpha=.5)+
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf), color='red', size=1)+
  #geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf), color='steelblue', size=1)+
  #ylim(0,.5)+
  ylab('BAF ProCESS')

dr_races = CNAqc:::blank_genome() + 
  geom_point(data = joint_table_snps_shifted %>% sample_n(size = n()/5), 
             aes(x=pos, y=DR), size = .1, alpha=.5)+
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_dr), color='red', size=1)+
  ylab('DR ProCESS')+
  ylim(0,4)

baf_comparison=CNAqc:::blank_genome() + 
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf+.02), color ='goldenrod', size=1)+
  geom_segment(data = inferred_baf_dr, aes(x=from, xend=to, y=inferred_baf-.02), color ='steelblue', size=1)+
  ylab('BAF') + 
  ggtitle('BAF: inferred (ascat, blue) vs simulated (ProCESS, yellow)')

dr_comparison=CNAqc:::blank_genome() + 
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_dr+.07), color ='goldenrod', size=1)+
  geom_segment(data = inferred_baf_dr, aes(x=from, xend=to, y=inferred_dr-.07), color ='steelblue', size=1)+
  ylab('DR') + 
  ylim(0,4) + 
  ggtitle('DR: inferred (ascat, blue) vs simulated (ProCESS, yellow)')

cna_calls_comparison = CNAqc:::blank_genome() + 
  geom_rect(data = joint_segmentation_shifted_longer,
            aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=is_match)) +
  geom_segment(data=joint_segmentation_shifted_longer %>% 
                 mutate(Value = case_when(
                   Type == 'TRUE_Major1' ~ Value + .09,
                   Type == 'TRUE_minor1' ~ Value - .09,
                   Type == 'INFERRED_Major1' ~ Value + .04,
                   Type == 'INFERRED_minor1' ~ Value - .04,
                   Type == 'TRUE_Major2' ~ Value + .09,
                   Type == 'TRUE_minor2' ~ Value - .09
                   )), 
               aes(x=from, xend=to, y=Value, color=Type), size=1)+
  scale_color_manual(values=color_by_state)+
  scale_fill_manual(values=fill_by_match)+
  labs(fill = "", color = "")+
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.text=element_text(size=12))


st = 
  'AAACCC
   BBBDDD
   EEEFFF
   GGGGGG
   GGGGGG'

report = patchwork::wrap_plots(
  baf_ascat, dr_ascat, baf_races, dr_races, baf_comparison, dr_comparison, cna_calls_comparison,
  design = st
) + patchwork::plot_annotation(
  title = element_text(paste0('CNA validation of sample: ', sample_id)),
  subtitle = element_text(paste0(
    'Proportion of genome inferred correctly: ', round(ascat_correctness,2)*100,'%',
    '\nTrue purity: ',purity_number,' Inferred purity: ', purity_ploidy$AberrantCellFraction,
    '\nTrue ploidy: ',round(ploidy,2),' Inferred ploidy: ', round(purity_ploidy$Ploidy,2),
    '\nNumber of subclonal segments: ', joint_segmentation_shifted %>% filter(!is.na(TRUE_Major1)) %>% 
      filter(!is.na(TRUE_Major2),TRUE_Major1!=TRUE_Major2,TRUE_minor1!=TRUE_minor2) %>% nrow())
)) 

outdir <- paste0(data_dir,spn_id,"/validation/cna/",spn_id,"/",coverage,"_",purity,'/',caller,"/", sample_id,"/")
dir.create(outdir, recursive = T, showWarnings = F)
ggsave(report, file = paste0(outdir,'report.png'), height = 12, width = 12)

reportdir <- paste0(data_dir,spn_id,"/validation/cna/report/")
dir.create(reportdir, recursive = T, showWarnings = F)
filename <- paste(spn_id,coverage,purity, caller, sample_id, sep='_')
file_path <- file.path(reportdir, filename)
ggsave(report, file = paste0(file_path,'.png'), height = 12, width = 12)

CNA_validation_summmary = list(
  'data' = joint_segmentation,
  'caller results' = CNA_ascat,
  'ProCESS results' = CNA_races,
  'SNPs subset' = joint_table_snps_shifted,
  'proportion of correctly inferred genome'= ascat_correctness,
  'true purity' = purity_number,
  'inferred purity' = purity_ploidy$AberrantCellFraction,
  'true ploidy' = ploidy,
  'inferred ploidy' = purity_ploidy$Ploidy,
  'simulated purity' = as.double(strsplit(purity, 'p')[[1]]),
  'state' = state
)

saveRDS(CNA_validation_summmary, file = paste0(outdir, 'metrics.rds'))
message("Report saved for combination: purity=", purity, ", cov=", coverage)




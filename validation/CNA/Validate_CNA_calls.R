options(bitmapType='cairo')
library(dplyr)
library(ProCESS)
library(optparse)
library(tidyr)
library(ggplot2)
library(future.apply)
library(progressr)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/sarek_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/antonelloa/ProCESS-examples/validation/CNA/utils.R")

## Installing sequenza
# dependencies ‘pbapply’, ‘squash’, ‘iotools’, ‘seqminer’, ‘copynumber’
# install.packages("pbapply")
# install.packages("squash")
# install.packages("iotools")
# install.packages("seqminer")
# install.packages("/orfeo/cephfs/scratch/cdslab/antonelloa/copynumber_1.34.0.tar", repos = NULL, type = "source")
# install.packages("/orfeo/cephfs/scratch/cdslab/antonelloa/sequenza_3.0.0.tar", repos = NULL, type = "source")

############ Parse command-line arguments
option_list <- list(make_option(c("--sample_id"), type = "character", default = 'SPN03_1.1'),
		                make_option(c("--spn_id"), type = "character", default = 'SPN03'),
                    make_option(c("--purity"), type = "character", default = '0.6'),
                    make_option(c("--coverage"), type = "character", default = '50'),
                    make_option(c("--purity_th"), type = "character", default = '.1'),
                    make_option(c("--correct_th"), type = "character", default = '.6')
		                # make_option(c("--caller"), type = "character", default = 'sequenza')
		                )
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'

sample_id = opt$sample_id
spn_id = opt$spn_id
coverage = paste0(opt$coverage)
purity = paste0(opt$purity)
purity_th = opt$purity_th
correct_th = opt$correct_th
# caller = opt$caller


#gender_file <- read.table(file = paste0(data_dir,spn_id,"/process/subject_gender.txt"),header = FALSE,col.names = "gender")
gender <- get_process_gender(spn = spn_id)

if (gender=="XX"){
  chromosomes = c(paste0('chr',1:22), 'chrX')
} else {
  chromosomes = c(paste0('chr',1:22), 'chrX', 'chrY')
}

############ Load data
#### ProCESS data

ProCESS_output = read_ProCESS(spn_id,sample_id,coverage,purity)
CNA_ProCESS = ProCESS_output[["CNA"]] 
snps = ProCESS_output[["snps"]] 

message("Reading ProCESS data")

#### ASCAT data

ASCAT_output = read_ASCAT(spn_id,sample_id,coverage,purity)
CNA_ascat = ASCAT_output[["CNA"]] 
BAF_DR_ascat = ASCAT_output[["BAF_DR"]] 
purity_ploidy_ascat = ASCAT_output[["purity_ploidy"]]

message("Reading ASCAT data")

#### Sequenza data

Sequenza_output = read_Sequenza(spn_id,sample_id,coverage,purity, seq_dir='/orfeo/cephfs/scratch/cdslab/antonelloa/ProCESS-examples/validation/CNA')
CNA_sequenza = Sequenza_output[["CNA"]] 
BAF_DR_sequenza = Sequenza_output[["BAF_DR"]] 
purity_ploidy_sequenza = Sequenza_output[["purity_ploidy"]] 
message("Reading Sequenza data")

#### CNVkit data

CNVkit_output = read_CNVkit(spn_id,sample_id,coverage,purity)
CNA_cnvkit = CNVkit_output[["CNA"]]
DR_cnvkit = CNVkit_output[["DR"]]

message("Reading CNVkit data")

############ Process data
message("Create joint table ProCESS and ASCAT calls") 
joint_segmentation_ascat = create_joint_segmentation(CNA_ProCESS, CNA_target=CNA_ascat, caller='ascat', chromosomes)
joint_segmentation_ascat_long = joint_segmentation_ascat[['joint_segmentation_long']]
message("Create joint table ProCESS and Sequenza calls") 
joint_segmentation_sequenza = create_joint_segmentation(CNA_ProCESS, CNA_target=CNA_sequenza, caller='sequenza', chromosomes)
joint_segmentation_sequenza_long = joint_segmentation_sequenza[['joint_segmentation_long']]
message("Create joint table ProCESS and CNVkit calls") 
joint_segmentation_cnvkit = create_joint_segmentation(CNA_ProCESS, CNA_target=CNA_cnvkit, caller='cnvkit', chromosomes)
joint_segmentation_cnvkit_long = joint_segmentation_cnvkit[['joint_segmentation_long']]

############ Compute CNA correctness
# Percentage of the genome correctly called (call is considered correct if the major and minor allele match, event if the match is reffered only to the most prevalent subclone)
message("Compute metrics")

ascat_correctness = compute_correctness(joint_segmentation_ascat_long) 
sequenza_correctness = compute_correctness(joint_segmentation_sequenza_long) 
cnvkit_correctness =  compute_correctness(joint_segmentation_cnvkit_long) 

purity_correctness_ascat = purity_ploidy_ascat$AberrantCellFraction - as.double(purity)
purity_correctness_sequenza = mean(purity_ploidy_sequenza$cellularity) - as.double(purity)

# Compute real ploidy
ploidy = compute_true_ploidy(CNA_ProCESS)


# Theoretical vs inferred BAF and DR
th_vs_inf_BAF_DR = theoretical_vs_inferred_BAF_DR(CNA_ProCESS, CNA_ascat, purity_number=as.double(purity),ploidy=ploidy,purity_ploidy=purity_ploidy_ascat)
simulated_baf_dr = th_vs_inf_BAF_DR[["simulated_baf_dr"]]
inferred_baf_dr_ascat = th_vs_inf_BAF_DR[["inferred_baf_dr_ascat"]]

## Segmentation correctness
# find out what proportion of each chromosome is not covered by the segmentation
# breakpoints:
#     - find matching BP and compute distance 
#     - find missed/addes BP
segmentation_ascat = segmentation_analysis(CNA_ProCESS, CNA_ascat, chromosomes, th=2e7)
segmentation_sequenza = segmentation_analysis(CNA_ProCESS, CNA_sequenza, chromosomes, th=2e7)
segmentation_cnvkit = segmentation_analysis(CNA_ProCESS, CNA_cnvkit, chromosomes, th=2e7)

breakpoints_ascat = segmentation_ascat[["breakpoints"]]
seg_summary_ascat = segmentation_ascat[["summary"]]
breakpoints_sequenza = segmentation_sequenza[["breakpoints"]]
seg_summary_sequenza = segmentation_sequenza[["summary"]]
breakpoints_cnvkit = segmentation_cnvkit[["breakpoints"]]
seg_summary_cnvkit = segmentation_cnvkit[["summary"]]
  

########### Plots
message("Generate plots")
color_by_state = c("TRUE_Major1"=alpha('firebrick', 1),
                   "TRUE_minor1"=alpha('#000080ff',1), 
                   "TRUE_Major2"=alpha('#ff00abb3'),
                   "TRUE_minor2"=alpha('slateblue', 1),
                   "INFERRED_Major1"=alpha('firebrick', .5),
                   "INFERRED_minor1"=alpha('#000080ff',.5),
                   'TRUE_CN'= alpha('firebrick', 1),
                   'INFERRED_CN'=alpha('#ff00abb3'))

fill_by_match = c('complete match'= alpha('gainsboro', .03),
                  'undetected subclone'= alpha('goldenrod', .08),
                  'close'= alpha('goldenrod', .08),
                  'no match' = alpha('indianred', .08)
)

fill_by_overlap <- colorRampPalette(c("firebrick", "goldenrod", "forestgreen"))(11)
names(fill_by_overlap) = as.character(seq(0, 100, by=10))

### Plots ASCAT
baf_ascat = CNAqc:::blank_genome() + 
  geom_point(data = BAF_DR_ascat , #%>% mutate(BAF = ifelse(BAF>.5, 1-BAF,BAF))
             aes(x=Position, y=BAF), size = .1, alpha=1)+
  geom_segment(data = inferred_baf_dr_ascat, aes(x=from, xend=to, y=inferred_baf), color='red', size=1)+
  #geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf), color='steelblue', size=1)+
  #ylim(0,.5)+
  ylab('BAF ascat')

dr_ascat = CNAqc:::blank_genome() + 
  geom_point(data = BAF_DR_ascat, 
             aes(x=Position, y=DR), size = .1, alpha=1)+
  geom_segment(data = inferred_baf_dr_ascat, aes(x=from, xend=to, y=inferred_dr), color='red', size=1)+
  ylab('DR ascat')+
  ylim(0,4)

baf_races = CNAqc:::blank_genome() + 
  geom_point(data = snps, #%>% mutate(BAF = ifelse(BAF>.5, 1-BAF,BAF)) 
             aes(x=pos, y=BAF), size = .1, alpha=.5)+
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf), color='red', size=1)+
  #geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf), color='steelblue', size=1)+
  #ylim(0,.5)+
  ylab('BAF ProCESS')

dr_races = CNAqc:::blank_genome() + 
  geom_point(data = snps, 
             aes(x=pos, y=DR), size = .1, alpha=.5)+
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_dr), color='red', size=1)+
  ylab('DR ProCESS')+
  ylim(0,4)

baf_comparison=CNAqc:::blank_genome() + 
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf+.02), color ='goldenrod', size=1)+
  geom_segment(data = inferred_baf_dr_ascat, aes(x=from, xend=to, y=inferred_baf-.02), color ='steelblue', size=1)+
  ylab('BAF') + 
  ggtitle('BAF: inferred (ascat, blue) vs simulated (ProCESS, yellow)')

dr_comparison=CNAqc:::blank_genome() + 
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_dr+.07), color ='goldenrod', size=1)+
  geom_segment(data = inferred_baf_dr_ascat, aes(x=from, xend=to, y=inferred_dr-.07), color ='steelblue', size=1)+
  ylab('DR') + 
  ylim(0,4) + 
  ggtitle('DR: inferred (ascat, blue) vs simulated (ProCESS, yellow)')

cna_calls_comparison = CNAqc:::blank_genome() + 
  geom_rect(data = joint_segmentation_ascat_long,
            aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=is_match)) +
  geom_segment(data=joint_segmentation_ascat_long %>% 
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


segmentation_comparison = my_blank_genome() +
  geom_segment(data=joint_segmentation_ascat_long %>% 
                 mutate(Value = case_when(
                   Type == 'TRUE_Major1' ~ Value + .12,
                   Type == 'TRUE_minor1' ~ Value - .12,
                   Type == 'INFERRED_Major1' ~ Value + .06,
                   Type == 'INFERRED_minor1' ~ Value - .06,
                   Type == 'TRUE_Major2' ~ Value + .12,
                   Type == 'TRUE_minor2' ~ Value - .12
                 )), 
               aes(x=from, xend=to, y=Value, color=Type), size=1)+
  scale_color_manual(values=color_by_state)+
  geom_vline(data = shift_segments(breakpoints_ascat %>% rename(from=og_coord, to=coord)) %>% rename(og_coord=from, coord=to) %>% filter(state=='matching'),
             aes(xintercept=og_coord), color='forestgreen', size=.3, linetype='dotted')+
  geom_vline(data = shift_segments(breakpoints_ascat %>% rename(from=og_coord, to=coord)) %>% rename(og_coord=from, coord=to) %>% filter(state=='missed'),
             aes(xintercept=og_coord), color='indianred', size=.3, linetype='dotted')+
  labs(fill = "", color = "")+
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.text=element_text(size=12))+theme(legend.position = 'none')+ggtitle('Breakpoints')
  


st = 
  'AAACCC
   BBBDDD
   EEEFFF
   GGGGGG
   GGGGGG
   HHHHHH'

report_ascat = patchwork::wrap_plots(
  baf_ascat, dr_ascat, baf_races, dr_races, baf_comparison, dr_comparison, cna_calls_comparison,segmentation_comparison,
  design = st
) + patchwork::plot_annotation(
  title = element_text(paste0('ASCAT CNA validation of sample: ', sample_id)),
  subtitle = element_text(paste0(
    'Proportion of genome inferred correctly: ', round(ascat_correctness,2)*100,'%',
    '\nTrue purity: ',purity,' Inferred purity: ', purity_ploidy_ascat$AberrantCellFraction,
    '\nTrue ploidy: ',round(ploidy,2),' Inferred ploidy: ', round(purity_ploidy_ascat$Ploidy,2),
    '\nNumber of subclonal segments: ', joint_segmentation_ascat[["joint_segmentation"]] %>% filter(!is.na(TRUE_Major1)) %>% 
      filter(!is.na(TRUE_Major2),TRUE_Major1!=TRUE_Major2,TRUE_minor1!=TRUE_minor2) %>% nrow(),
    '\nAverage breakpoint distance (matching inferred vs real): ',round(seg_summary_ascat$av_distance)
))) 

### Plots Sequenza
baf_sequenza = CNAqc:::blank_genome() + 
  geom_point(data = BAF_DR_sequenza %>% sample_n(size = n()/2) %>% filter(depth.tumor > 60), #%>% mutate(BAF = ifelse(BAF>.5, 1-BAF,BAF))
             aes(x=shifted_position, y=Bf), size = .1, alpha=1)+
  geom_segment(data = shift_segments(CNA_sequenza), aes(x=from, xend=to, y=BAF), color='red', size=1)+
  #geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf), color='steelblue', size=1)+
  #ylim(0,.5)+
  ylab('BAF sequenza')

dr_sequenza = CNAqc:::blank_genome() + 
  geom_point(data = BAF_DR_sequenza %>% sample_n(size = n()/2) %>% filter(depth.tumor > 60) %>% rowwise() %>% 
               mutate(normalised_dr=exp(log(depth.ratio*(depth.normal/depth.tumor)))), 
             aes(x=shifted_position, y=normalised_dr), size = .1, alpha=1)+
  geom_segment(data = shift_segments(CNA_sequenza), aes(x=from, xend=to, y=DR), color='red', size=1)+
  ylab('DR sequenza')+
  ylim(0,4)

baf_comparison_sequenza=CNAqc:::blank_genome() + 
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_baf+.02), color ='goldenrod', size=1)+
  geom_segment(data = shift_segments(CNA_sequenza), aes(x=from, xend=to, y=BAF-.02), color ='steelblue', size=1)+
  ylab('BAF') + 
  ggtitle('BAF: inferred (sequenza, blue) vs simulated (ProCESS, yellow)')

dr_comparison_sequenza=CNAqc:::blank_genome() + 
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_dr+.07), color ='goldenrod', size=1)+
  geom_segment(data = shift_segments(CNA_sequenza), 
               aes(x=from, xend=to, y=DR-.07), color ='steelblue', size=1)+
  ylab('DR') + 
  ylim(0,4) + 
  ggtitle('DR: inferred (sequenza, blue) vs simulated (ProCESS, yellow)')

cna_calls_comparison_sequenza = CNAqc:::blank_genome() + 
  geom_rect(data = joint_segmentation_sequenza_long,
            aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=is_match)) +
  geom_segment(data=joint_segmentation_sequenza_long %>% 
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

segmentation_comparison_sequenza = my_blank_genome() +
  geom_segment(data=joint_segmentation_sequenza_long %>% 
                 mutate(Value = case_when(
                   Type == 'TRUE_Major1' ~ Value + .12,
                   Type == 'TRUE_minor1' ~ Value - .12,
                   Type == 'INFERRED_Major1' ~ Value + .06,
                   Type == 'INFERRED_minor1' ~ Value - .06,
                   Type == 'TRUE_Major2' ~ Value + .12,
                   Type == 'TRUE_minor2' ~ Value - .12
                 )), 
               aes(x=from, xend=to, y=Value, color=Type), size=1)+
  scale_color_manual(values=color_by_state)+
  geom_vline(data = shift_segments(breakpoints_sequenza %>% rename(from=og_coord, to=coord)) %>% rename(og_coord=from, coord=to) %>% filter(state=='matching'),
             aes(xintercept=og_coord), color='forestgreen', size=.3, linetype='dotted')+
  geom_vline(data = shift_segments(breakpoints_sequenza %>% rename(from=og_coord, to=coord)) %>% rename(og_coord=from, coord=to) %>% filter(state=='missed'),
             aes(xintercept=og_coord), color='indianred', size=.3, linetype='dotted')+
  labs(fill = "", color = "")+
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.text=element_text(size=12))+theme(legend.position = 'none')+ggtitle('Breakpoints')


st = 
  'AAACCC
   BBBDDD
   EEEFFF
   GGGGGG
   GGGGGG
   HHHHHH'

report_sequenza = patchwork::wrap_plots(
  baf_sequenza, dr_sequenza, baf_races, dr_races, baf_comparison_sequenza, dr_comparison_sequenza, cna_calls_comparison_sequenza,segmentation_comparison_sequenza,
  design = st
) + patchwork::plot_annotation(
  title = element_text(paste0('Sequenza CNA validation of sample: ', sample_id)),
  subtitle = element_text(paste0(
    'Proportion of genome inferred correctly: ', round(sequenza_correctness,2)*100,'%',
    '\nTrue purity: ',purity,' Inferred purity: ', round(mean(purity_ploidy_sequenza$cellularity),3),
    '\nTrue ploidy: ',round(ploidy,2),' Inferred ploidy: ', mean(round(purity_ploidy_sequenza$ploidy.estimate,2)),
    '\nNumber of subclonal segments: ', joint_segmentation_sequenza[["joint_segmentation"]] %>% filter(!is.na(TRUE_Major1)) %>% 
      filter(!is.na(TRUE_Major2),TRUE_Major1!=TRUE_Major2,TRUE_minor1!=TRUE_minor2) %>% nrow(),
    '\nAverage breakpoint distance (matching inferred vs real): ',round(seg_summary_ascat$av_distance))
  )) 


### Plots CNVkit
dr_cnvkit = CNAqc:::blank_genome() + 
  geom_segment(data = DR_cnvkit %>% sample_n(size = n()/2), 
             aes(x=start,xend=end, y=weight, yend=weight), size = 2, alpha=1)+
  geom_segment(data = shift_segments(CNA_cnvkit), aes(x=from, xend=to, y=CN/2), color='red', size=1)+
  ylab('DR CNVkit')+
  ylim(0,4)

#dr_races 

dr_comparison_cnvkit=CNAqc:::blank_genome() + 
  geom_segment(data = simulated_baf_dr, aes(x=from, xend=to, y=simulated_dr+.07), color ='goldenrod', size=1)+
  geom_segment(data = shift_segments(CNA_cnvkit), aes(x=from, xend=to, y=CN/2-.07), color ='steelblue', size=1)+
  ylab('DR') + 
  ylim(0,4) + 
  ggtitle('DR: inferred (cnvkit, blue) vs simulated (ProCESS, yellow)')

cna_calls_comparison_cnvkit = CNAqc:::blank_genome() + 
  geom_rect(data = joint_segmentation_cnvkit_long,
            aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=is_match)) +
  geom_segment(data=joint_segmentation_cnvkit_long %>% 
                 mutate(Value = case_when(
                   Type == 'TRUE_CN' ~ Value + .09,
                   Type == 'INFERRED_CN' ~ Value - .09
                 )), 
               aes(x=from, xend=to, y=Value, color=Type), size=1)+
  scale_color_manual(values=color_by_state)+
  scale_fill_manual(values=fill_by_match)+
  labs(fill = "", color = "")+
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.text=element_text(size=12))+
  ylab('Total CN')

segmentation_comparison_cnvkit = my_blank_genome() +
  geom_segment(data=joint_segmentation_sequenza_long %>% 
                 mutate(Value = case_when(
                   Type == 'TRUE_Major1' ~ Value + .12,
                   Type == 'TRUE_minor1' ~ Value - .12,
                   Type == 'INFERRED_Major1' ~ Value + .06,
                   Type == 'INFERRED_minor1' ~ Value - .06,
                   Type == 'TRUE_Major2' ~ Value + .12,
                   Type == 'TRUE_minor2' ~ Value - .12
                 )), 
               aes(x=from, xend=to, y=Value, color=Type), size=1)+
  scale_color_manual(values=color_by_state)+
  geom_vline(data = shift_segments(breakpoints_cnvkit %>% rename(from=og_coord, to=coord)) %>% rename(og_coord=from, coord=to) %>% filter(state=='matching'),
             aes(xintercept=og_coord), color='forestgreen', size=.3, linetype='dotted')+
  geom_vline(data = shift_segments(breakpoints_cnvkit %>% rename(from=og_coord, to=coord)) %>% rename(og_coord=from, coord=to) %>% filter(state=='missed'),
             aes(xintercept=og_coord), color='indianred', size=.3, linetype='dotted')+
  labs(fill = "", color = "")+
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.text=element_text(size=12))+theme(legend.position = 'none')+ggtitle('Breakpoints')


st = 
  'AAABBB
   CCC###
   DDDDDD
   DDDDDD
   EEEEEE'

report_cnvkit = patchwork::wrap_plots(
  dr_cnvkit, dr_races, dr_comparison_cnvkit, cna_calls_comparison_cnvkit,segmentation_comparison_cnvkit,
  design = st
) + patchwork::plot_annotation(
  title = element_text(paste0('CNVkit CNA validation of sample: ', sample_id)),
  subtitle = element_text(paste0(
    'Proportion of genome inferred correctly: ', round(cnvkit_correctness,2)*100,'%',
    '\nAverage breakpoint distance: ',round(seg_summary_cnvkit$av_distance)
  ))) 



### Save reports 
outdir <- paste0(data_dir,spn_id,"/validation/cna/",spn_id,"/",coverage,"_",purity,'/',sample_id,'/')
dir.create(outdir, recursive = T, showWarnings = F)
# ascat
ggsave(report_ascat, file = paste0(outdir,'ascat_report.png'), height = 12, width = 12)
# cnvkit
ggsave(report_cnvkit, file = paste0(outdir,'cnvkit_report.png'), height = 7, width = 12)

# reportdir <- paste0(data_dir,spn_id,"/validation/cna/report/")
# dir.create(reportdir, recursive = T, showWarnings = F)
# filename <- paste(spn_id,coverage,purity, caller, sample_id, sep='_')
# file_path <- file.path(reportdir, filename)
# ggsave(report, file = paste0(file_path,'.png'), height = 12, width = 12)

# if (purity_th < abs(purity_correctness_ascat) & 
#     purity_th < abs(purity_correctness_sequenza) &
#     correct_th < ascat_correctness & 
#     correct_th < cnvkit_correctness &
#     ){state = 'PASS'}else{state='FAIL'}

CNA_validation_summmary = list(
  #'data' = joint_segmentation,
  'ascat results' = CNA_ascat,
  'cnvkit results' = CNA_cnvkit,
  'ProCESS results' = CNA_races,
  'SNPs subset' = joint_table_snps_shifted,
  'proportion of correctly inferred genome (ascat)'= ascat_correctness,
  'proportion of correctly inferred genome (cnvkit)'= cnvkit_correctness,
  'true purity' = purity_number,
  'inferred purity (ascat)' = purity_ploidy$AberrantCellFraction,
  'true ploidy' = ploidy,
  'inferred ploidy (ascat)' = purity_ploidy$Ploidy,
  'simulated purity' = as.double(strsplit(purity, 'p')[[1]]),
  'average bp distance (ascat)'= mean(filtered_best_segments$distance),
  'average bp distance (cnvkit)'= mean(filtered_best_segments$distance),
  'state' = state
)
#outdir <- paste0(data_dir,spn_id,"/validation/cna/",spn_id,"/",coverage,"_",purity,'/sample_id/')
saveRDS(CNA_validation_summmary, file = paste0(outdir, 'metrics.rds'))
message("Report saved for combination: purity=", purity, ", cov=", coverage)




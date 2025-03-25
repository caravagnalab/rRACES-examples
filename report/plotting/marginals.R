library(dplyr)
library(rRACES)
library(ggplot2)
seq_res <- readRDS("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/races/purity_0.9/seq_results_muts_merged_coverage_100x.rds")
phylo_forest <- load_phylogenetic_forest("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/races/phylo_forest.sff")
sample_forest <- load_samples_forest("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/races/sample_forest.sff")

# driver_chroms <- phylo_forest$get_driver_mutations() %>% 
#   filter(type=="SID") %>% 
#   pull(chr)
# driver_table <-annotate_drivers(phylo_forest) %>% 
#   filter(type=="SID") %>% 
#   rename(from=pos)
data <- seq_res %>%
  filter(classes!="germinal") %>% 
  # filter(chr %in% driver_chroms) %>%
  filter(!stringr::str_detect(causes, 'errors')) %>%
  rRACES::seq_to_long()

# Set desired fraction of total mutations to sample
fraction_to_sample <- 0.1  # 30%

# Compute number of samples per chromosome while keeping proportions
sample_sizes <- data %>%
  ungroup() %>% 
  count(chr) %>%
  mutate(sample_size = round(n * fraction_to_sample))

# Sample mutations proportionally from each chromosome
df_sampled <- data %>%
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



mutations = df_sampled %>% as_tibble()%>%
  mutate(mutation_id=paste(chr, from, to, ref, alt, sep="_"), chr_pos=from)%>%
  select(mutation_id, everything())


mutations = mutations %>%
  filter(rowSums(select(., ends_with(".NV")) != 0) > 0)

model_df = mutations %>%
  filter(!stringr::str_detect(causes, 'errors')) %>%
  select(mutation_id, contains("VAF"))


mutation_ids <- unique(model_df$mutation_id)

labels <- c()
tot_samples <- 3
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

mutations <- as.data.frame(mutations)
mutations_with_cell = mutations %>%
  filter(classes!="germinal") %>% 
  filter(!stringr::str_detect(causes, 'errors')) %>% 
  rowwise() %>%
  mutate(cell_id=phylo_forest$get_first_occurrences(Mutation(
    chr, chr_pos, ref, alt
  ))[[1]]) %>%
  ungroup()

cells_labels = mutations_with_cell %>% 
  select(mutation_id, cell_id) %>%
  left_join(model_df) %>% 
  group_by(cell_id) %>% 
  summarise(label_list=list(label)) %>% 
  rowwise() %>% 
  mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
  ungroup() %>% 
  select(-label_list)

final_labels = sample_forest$get_nodes() %>% as_tibble() %>% 
  left_join(cells_labels)
palette_labels <- RColorBrewer::brewer.pal(n = length(unique(labels)),name = "Set2")
names(palette_labels) <-unique(labels)
pl_sticks <- plot_sticks(sample_forest, labels=final_labels,cls = palette_labels) %>%
  annotate_forest(sample_forest, samples=TRUE, drivers=TRUE)

pl = patchwork::wrap_plots(pl_sticks, design="AAA\nAAA")





plot_vaf_marginal <- lapply(pairwise_comb,function(x){
  s1 <- x[1]
  s2 <- x[2]
  p = model_df %>% 
    ggplot(aes(x =.data[[vaf[s1]]], y = .data[[vaf[s2]]],color=(label)))+
    geom_point()+
    # ggrepel::geom_label_repel(aes(label = driver_label),
    #                           na.rm = TRUE, box.padding = 0.5,
    #                           color = "black",size=2
    # ) + 
    xlab(s1) +
    ylab(s2)+
    scale_color_manual(values = palette_labels)+
    CNAqc:::my_ggplot_theme()+
    theme(legend.position = "none")
  
  return(p)
})
pl_mrg <- wrap_plots(plot_vaf_marginal,ncol = 1)

# p_my_lab <- patchwork::wrap_plots(plot_vaf_marginal, design="AAAB\nAAAC\nAAAD")
# p_my_lab
patchwork::wrap_plots(list(pl,pl_mrg), design="AAAB\nAAAB\nAAAB")

patchwork::wrap_plots(list(pl,pl_mrg,pl_hist), design="AAAB\nAAAB\nAAAB\nCCCC") & theme(legend.position = "bottom")
####################################################################



drivers_long <- s_seq_long %>% filter(classes=="driver") %>% 
  mutate(mutation_id=paste(chr,from,ref,alt,sep=":"))

ggplot(drivers_long) +
  geom_line(aes(x = sample_name, y = VAF, group = mutation_id, color = causes)) +
  geom_point(aes(x = sample_name, y = VAF, color = causes)) +
  labs(y = "VAF") +
  theme_bw() +
  scale_color_manual(values=get_clone_map(sample_forest))+
  ylim(0,1)





####################################################################
km = kmeans(model_df %>% select(contains("VAF")), centers=8)
model_df$label_kmeans = km$cluster
color_palette = c(
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#a65628",
  "#FFD700",  "#999999", "#000000", "#f781bf", # First 10 colors (Set1)
  "#46f0f0", "#f032e6", "#bcf60c", "#fabed4", "#008080", "#e6beff",
  "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",
  "#000075", "#808080", "#d3a6f3", "#ff9cdd", "#73d7b0"
)
cells_labels = mutations_with_cell %>% 
  select(mutation_id, cell_id) %>%
  left_join(model_df) %>% 
  group_by(cell_id) %>% 
  summarise(label_list=list(label_kmeans)) %>% 
  rowwise() %>% 
  mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
  ungroup() %>% 
  select(-label_list)

final_labels = sample_forest$get_nodes() %>% as_tibble() %>% 
  left_join(cells_labels)



pl_sticks <- plot_sticks(sample_forest, labels=final_labels, cls=color_palette) %>%
  annotate_forest(sample_forest, samples=TRUE, drivers=TRUE)

pl = patchwork::wrap_plots(pl_sticks, design="AAA\nAAA")
pl

p1 <-model_df %>% 
  ggplot(aes(x=SPN01_1.1.VAF,y=SPN01_1.2.VAF,color=factor(label_kmeans)))+
  geom_point()+
  ggrepel::geom_label_repel(aes(label = driver_label,colour = "black"),
                            na.rm = TRUE, box.padding = 0.5
                            # nudge_y = c(-0.2)
  ) + 
  scale_color_manual(values = color_palette)+
  CNAqc:::my_ggplot_theme()+
  theme(legend.position = "none")
p1_marg <- ggMarginal(p1, type = "density", groupColour = TRUE, groupFill = TRUE)

p2 <- model_df %>% 
  ggplot(aes(x=SPN01_1.2.VAF,y=SPN01_1.3.VAF,color=factor(label_kmeans)))+
  geom_point()+
  ggrepel::geom_label_repel(aes(label = driver_label,colour = "black"),
                            na.rm = TRUE, box.padding = 0.5
                            # nudge_y = c(-0.2)
  ) + 
  scale_color_manual(values = color_palette)+
  CNAqc:::my_ggplot_theme()+
  theme(legend.position = "none")
p2_marg <- ggMarginal(p2, type = "density", groupColour = TRUE, groupFill = TRUE)

p3 <- model_df %>% 
  ggplot(aes(x=SPN01_1.1.VAF,y=SPN01_1.3.VAF,color=factor(label_kmeans)))+
  geom_point()+
  ggrepel::geom_label_repel(aes(label = driver_label,colour = "black"),
                            na.rm = TRUE, box.padding = 0.5
                            # nudge_y = c(-0.2)
  ) + 
  scale_color_manual(values = color_palette)+
  CNAqc:::my_ggplot_theme()+
  theme(legend.position = "none")
p3_marg <- ggMarginal(p3, type = "density", groupColour = TRUE, groupFill = TRUE)
p_kmeans <- patchwork::wrap_plots(list(pl,p1_marg,p2_marg,p3_marg), design="AAAB\nAAAC\nAAAD")
p_kmeans



library(tidyverse)
library(ggalluvial)

# Step 1: Classify mutations
sankey_data <- model_df %>%
  select(c("SPN01_1.3.VAF","SPN01_1.1.VAF")) %>% 
  filter(rowSums(select(., ends_with(".VAF")) != 0) > 0) %>% 
  mutate(Category = case_when(
    SPN01_1.1.VAF > 0 & SPN01_1.3.VAF > 0 ~ "Clonal Shared",
    SPN01_1.1.VAF > 0 & SPN01_1.3.VAF == 0 ~ "Private Primary",
    SPN01_1.1.VAF == 0 & SPN01_1.3.VAF > 0 ~ "Private Relapse",
    TRUE ~ "Missing"
  )) %>%
  count(Category) %>%
  pivot_longer(cols = n, names_to = "Count", values_to = "Size")

# Step 2: Prepare for ggalluvial
sankey_long <- sankey_data %>%
  mutate(Primary = recode(Category,
                          "Clonal Shared" = "Clonal Primary",
                          "Private Primary" = "Tail Primary",
                          "Private Relapse" = "Missing",
                          "Missing" = "Missing"),
         Relapse = recode(Category,
                          "Clonal Shared" = "Clonal Relapse",
                          "Private Primary" = "Missing",
                          "Private Relapse" = "Tail Relapse",
                          "Missing" = "Missing")) %>%
  select(Primary, Relapse, Size)

# Step 3: Create Sankey Diagram
ggplot(sankey_long, aes(axis1 = Primary, axis2 = Relapse, y = Size)) +
  geom_alluvium(aes(fill = Primary), width = 0.3, alpha = 0.7) +
  geom_stratum(width = 0.3, fill = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5) +
  theme_minimal() +
  labs(title = "Sankey Diagram of Mutation Categories",
       x = "Primary",
       y = "Group Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


purity <- params$sequencing$purity
table <- samples_table(snapshot=params$files$sim,
                       forest=params$files$sample_forest)

table <- table %>% 
  mutate(Normal_Cells = round(Total_Cells*purity,0)) %>% 
  mutate(All_Cells = Total_Cells+Normal_Cells) %>% 
  select("Sample_ID","Clone 3","Clone 4","Normal_Cells","All_Cells")

data_long <- table %>%
  pivot_longer(cols = c("Clone 3", "Clone 4", "Normal_Cells"),
               names_to = "Cell_Type",
               values_to = "Count")

# Compute proportions
data_long <- data_long %>%
  group_by(Sample_ID) %>%
  mutate(Prop = Count / sum(Count)) %>%
  ungroup()

# Create inner layer (Sample_ID) proportional segments
inner_data <- table %>%
  mutate(Prop = All_Cells / sum(All_Cells)) %>%
  select(Sample_ID, Prop) %>%
  mutate(Level = "Sample")

# Create outer layer (Cell types) with correct proportions
outer_data <- data_long %>%
  mutate(Level = "Cell_Type")

# Combine datasets
plot_data <- bind_rows(
  inner_data %>% rename(Category = Sample_ID),
  outer_data %>% rename(Category = Cell_Type)
)

# Assign levels to control donut rings
plot_data <- plot_data %>%
  mutate(Ring = ifelse(Level == "Sample", 1, 2))

# Compute segment positions
plot_data <- plot_data %>%
  arrange(Level, Category) %>%
  mutate(
    xmin = c(0, head(cumsum(Prop), -1)),
    xmax = cumsum(Prop),
    xmid = (xmin + xmax) / 2,  # Midpoint for text positioning
    ymid = ifelse(Ring == 1, 1, 2),  # Inner ring (1), Outer ring (2)
    ymin = ymid - 0.4, ymax = ymid + 0.4
  )

# Create the nested donut chart
ggplot(plot_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Category)) +
  geom_rect(color = "white") +
  geom_text(aes(x = xmid, y = ymid, label = Category), size = 5) +  # Fix text placement
  coord_polar(theta = "x") +
  theme_void() +
  theme(legend.position = "none")

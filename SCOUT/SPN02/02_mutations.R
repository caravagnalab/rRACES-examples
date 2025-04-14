# adding mutations

# Mutation generation ####
library(ProCESS)
library(dplyr)
library(ggplot2)

clone3_born = readRDS("/orfeo/cephfs/scratch/area/vgazziero/CDSlab/SPN02/results/clone3_clock.rds")
sampled_phylogeny = load_samples_forest("/orfeo/cephfs/scratch/area/vgazziero/CDSlab/SPN02/results/samples_forest.sff")
curr_dir = getwd()
setwd("/orfeo/cephfs/scratch/cdslab/shared/races/")

# setting the germline subject as columbian woman and tumor type to colonrectal cancer
m_engine <- MutationEngine(setup_code = "GRCh38", 
                context_sampling = 50, 
                tumour_type = "COAD",
			   tumour_study = "US", germline_subject = "HG01113")


m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = 1e-8, indel = 1e-9, CNA = 1e-11),
                    drivers = list("BRAF V600E"))

m_engine$add_mutant(mutant_name = "Clone 2", 
                    passenger_rates = c(SNV = 1e-8, indel = 1e-9, CNA = 1e-11),
                    drivers = list("PIK3CA E545K"))

m_engine$add_mutant(mutant_name = "Clone 3", 
                    passenger_rates = c(SNV = 1e-7, indel = 1e-8, CNA = 1e-10), # modelled as an hypermutant, so with higher passenger_rates
                    drivers = list(SNV("2", 47799065, "A")))

# indels and sbs have different coefficients summing (up to 1 for IDs and up to 1 for SBSs)
m_engine$add_exposure(
    time = 0, 
    coefficients = c(
        ID1 = 1,
        SBS1 = 0.5,
        SBS5 = 0.5)
)

m_engine$add_exposure(
    time = round(clone3_born),
    coefficients = c(
        SBS1 = 0.2, SBS5 = 0.4, SBS6 = 0.4, 
        ID1 = 0.5, ID7 = 0.5)
)

# Sys.time()
# add mutations on the forest with 1000 pre-neoplastic mutations
phylo_forest <- m_engine$place_mutations(sampled_phylogeny, 
    num_of_preneoplatic_SNVs = 800, 
    num_of_preneoplatic_indels = 200
)
Sys.time()
print("Mutations placed")
# save the phylogenetic forest in the file "phylo_forest.sff"
phylo_forest$save(paste(curr_dir, "phylo_forest.sff", sep = "/"))

# plotting
tree_plot = plot_forest(sampled_phylogeny) 
mut_forest = annotate_forest(tree_plot,
                phylo_forest,
                samples = TRUE,
                MRCAs = TRUE,
                exposures = T,
                facet_signatures = F,
                drivers = T,
                add_driver_label = T) +
            ggplot2::facet_wrap( ~ signature, scales = "free")

ggsave(filename = "/orfeo/cephfs/scratch/area/vgazziero/CDSlab/ProCESS/ProCESS-examples/SPN02/phylogenetic_forest_cycles.png", plot = mut_forest, width = 20, height = 15, bg = "white")
print("done everything")

exposure_time = plot_exposure_timeline(
  phylo_forest,
  linewidth = 0.8,
  emphatize_switches = TRUE
) 
ggsave(filename = "/orfeo/scratch/area/vgazziero/CDSlab/ProCESS/ProCESS-examples/SPN02/exposure_plot.png", plot = exposure_time, bg = "white")

all = mut_forest / exposure_time
ggsave(filename = "/orfeo/scratch/area/vgazziero/CDSlab/ProCESS/ProCESS-examples/SPN02/recap_phylo.png", plot = all)
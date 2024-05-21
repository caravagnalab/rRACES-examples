# adding mutations

# Mutation generation ####
library(rRACES)
library(dplyr)

clone3_born = readRDS("test/SPN02/clone3_clock.rds")
sampled_phylogeny = load_samples_forest("test/SPN02/samples_forest.sff")
curr_dir = getwd()
setwd("/orfeo/cephfs/scratch/cdslab/shared/races")
m_engine <- build_mutation_engine(setup_code = "GRCh38", 
                context_sampling = 50)
setwd(curr_dir)
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-11),
                    drivers = list("BRAF V600E"))

m_engine$add_mutant(mutant_name = "Clone 2", 
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-11),
                    drivers = list("PIK3CA E545K"))

m_engine$add_mutant(mutant_name = "Clone 3", 
                    passenger_rates = c(SNV = 1e-7, CNA = 1e-11), # modelled as an hypermutant, so with higher passenger_rates
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
    time = clone3_born,
    coefficients = c(
        SBS1 = 0.2, SBS5 = 0.4, SBS6 = 0.4, 
        ID1 = 0.5, ID7 = 0.5)
)

Sys.time()
# add mutations on the forest with 1000 pre-neoplastic mutations
phylo_forest <- m_engine$place_mutations(sampled_phylogeny, num_of_preneoplatic_mutations = 1000)
Sys.time()
print("Mutations placed")
# save the phylogenetic forest in the file "phylo_forest.sff"
phylo_forest$save("test/SPN02/phylo_forest.sff")

# plotting
tree_plot = plot_forest(sampled_phylogeny) 
mut_forest = annotate_forest(tree_plot,
                phylo_forest,
                samples = TRUE,
                MRCAs = TRUE,
                exposures = T,
                facet_signatures = T,
                drivers = T,
                add_driver_label = T) +
                ylim(0,10)

ggsave(filename = "test/SPN02/phylogenetic_forest_cycles.pdf", plot = mut_forest, width = 15, height = 10)
print("done everything")

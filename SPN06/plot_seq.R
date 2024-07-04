
# plot seq

#-------------------------------------------------------------------------------
#----------------------------------- plots -------------------------------------
#-------------------------------------------------------------------------------


# Mutation simulation
p1 <- plot_forest(forest = 0, highlight_sample = NULL)
p2 <- plot_sticks(forest = 0, labels = 0, cls = NULL)
p3 <- plot_exposure_timeline(phylogenetic_forest = 0, linewidth = 0.8, emphatize_switches = FALSE)

# Sequencing
p4 <- plot_BAF(seq_res = 0, sample = 0, chromosomes = NULL, cuts = c(0,1), N = 5000)
p5 <- plot_DR(seq_res = 0, sample = 0, chromosomes = NULL, N = 5000)
p6 <- plot_VAF(seq_res = 0, sample = 0, chromosomes = NULL, cuts = c(0,1), N = 5000)
p7 <- plot_VAF_histogram(seq_res = 0, chromosomes = NULL, samples = NULL, labels = NULL, cuts = c(0,1))
p8 <- plot_VAF_marginals(seq_res = 0, chromosomes = NULL, samples = NULL, labels = NULL, cuts = c(0,1))


plot_BAF(seq_res = seq_results, chromosomes = NULL, cuts = c(0,1), N = 5000)


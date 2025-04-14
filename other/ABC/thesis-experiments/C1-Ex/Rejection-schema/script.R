# Import librarys
library(EasyABC)
library(ProCESS)
library(ggplot2)

# Setting working directory path
path = getwd()
print(path)

# setting seed to make the experiment reproducible
set.seed(1)


############################################################################
# Creation of the directory where the plots will be saved
sub_dir <- "plots"

# check if sub directory exists 
if (file.exists(sub_dir)){
    
    # specifying the working directory
    setwd(file.path(path, sub_dir))
} else {
    
    # create a new sub directory inside
    # the main path
    dir.create(file.path(path, sub_dir))
    
    # specifying the working directory
    setwd(file.path(path, sub_dir))
}


############################################################################
# Simulation 1

# Setting parameters to generate data
# growth rate
s1.gr = 0.08567
# death rate
s1.dr = 0
# time units of evolution
s1.t=150
# Setting parameter for the ABC simulation
# number of simulations
nb_simul=6000
# tolerance rate (percentage of simulation accepted)
tr=0.05


############################################################################
# Data generation
# Create simulation
sim <- new(Simulation, 'sim')
sim$history_delta = 1

# Add the mutant A
sim$add_mutant('A', s1.gr, s1.dr)

# Place cell on the tissue
dim <- sim$get_tissue_size()
sim$place_cell('A', dim[1]/2, dim[2]/2)

# Run the simulation for t time units
sim$run_up_to_time(s1.t)

# Get the summary statistic of the simulation (number of cells)
s1.sum_stat_obs = c(sim$get_counts()$counts[1])


############################################################################
# Saving the Tissue and Muller plots
p=plot_tissue(sim)
ggsave("tissue_plot_sim1.png", plot=p, path = paste0(path, "/plots"))

p=plot_muller(sim)
ggsave("muller_plot_sim1.png", plot=p, path = paste0(path, "/plots"))


############################################################################
# Creation of the ProCESS model which reproduces the previous phenomenon and 
# which will be used for the ABC analysis

# It simply reproduces the previous steps
create_model <- function(t){
    model=function(x){
        if (x[1] > 0){
            tmp <- new(Simulation, 'tmp')
            tmp$add_mutant('A', x[1], 0)
            
            dim <- tmp$get_tissue_size()
            tmp$place_cell("A", dim[1]/2, dim[2]/2)
            tmp$run_up_to_time(tmp$get_clock() + t)
            
            c(tmp$get_counts()$counts[1])
        } else {
            c(0)
        }
    }
    print("The model has been created")
    return(model)
}
print("Modello creato")

# Instantiation of the model with the set t
s1.model = create_model(s1.t)

# We assume that we have no information on the growth rate and use a non-informative prior,
# therefore a uniform distribution
s1.prior = list(c("unif",0,0.1))

# Begin the ABC simulation
print("Inizio simulazione ABC 1")
s1.ABC_rej <- ABC_rejection(
    model=s1.model,
    prior=s1.prior,
    nb_simul=nb_simul,
    summary_stat_target=s1.sum_stat_obs,
    tol=tr)

print("Fine simulazione ABC 1")

############################################################################
# Prepare the output of the ABC simulation
s1.ABC_param = data.frame("p1" = s1.ABC_rej$param)

# Save the accepted value
p=ggplot(s1.ABC_param, aes(x=p1)) +
    
    geom_histogram(color="darkblue", fill="lightblue") +
    
    geom_vline(aes(xintercept = mean(p1)), color = "red", linetype="longdash") +
    
    geom_vline(aes(xintercept = mean(s1.gr)), color = "blue", linetype="longdash") +
    
    labs(y = "counts", x = expression(theta)) +
    
    theme_minimal()

ggsave("estimated_parameter_sim1.png", plot=p, path = paste0(path, "/plots"))

# Save tha accetped value in csv file
s1.tab <- cbind(s1.ABC_rej$param[,1])
s1.tab <- as.table(s1.tab)
write.table(s1.tab, file = "sim1_result.csv", sep = ",", col.names = NA,
            qmethod = "double")


print("Fine simulazione 1")

#######################################################################################################

############################################################################
# Simulation 2
print("Inizio simulazione 2")

# Set seed to make the experiment reproducible
set.seed(1)

# Setting parameters to generate data
# growth rate
s2.gr = 0.07
# death rate
s2.dr = 0.015
# time units of evloution
s2.t=150

# Setting parameter for the ABC simulation
# number of simulations
nb_simul=6000
# tolerance rate (percentage of simulation accepted)
tr=0.05


##########################################################################
# Data generation
# Create simulation
sim <- new(Simulation, 'sim')
sim$history_delta = 1

# Add the mutant A
sim$add_mutant('A', s2.gr, s2.dr)

# Place cell on the tissue
dim <- sim$get_tissue_size()
sim$place_cell('A', dim[1]/2, dim[2]/2)

# Run the simulation for other t time units
sim$run_up_to_time(s2.t)

# Get the summary statistic of the simulation (number of cells)
s2.sum_stat_obs = c(sim$get_counts()$counts[1])


############################################################################
# Saving the Tissue and Muller plots
p=plot_tissue(sim)
ggsave("tissue_plot_sim2.png", plot=p, path = paste0(path, "/plots"))

p=plot_muller(sim)
ggsave("muller_plot_sim2.png", plot=p, path = paste0(path, "/plots"))

############################################################################
# Creation of the ProCESS model which reproduces the previous phenomenon and 
# which will be used for the ABC analysis

# It simply reproduces the previous steps
create_model <- function(t){
    model=function(x){
        if (all(x > 0)){
            tmp <- new(Simulation, 'tmp')
            tmp$add_mutant('A', x[1], x[2])
            
            dim <- tmp$get_tissue_size()
            tmp$place_cell("A", dim[1]/2, dim[2]/2)
            tmp$run_up_to_time(tmp$get_clock() + t)
            
            c(tmp$get_counts()$counts[1])
        } else {
            c(0)
        }
    }
    print("The model has been created")
    return(model)
}

# Instantiation of the model with the set t
s2.model = create_model(s2.t)

# We assume that we have no information on the growth rate and use a non-informative prior,
# therefore a uniform distribution
s2.prior = list(c("unif",0.03,0.10), c("unif",0,0.06))

print("Inizio simulazione ABC 2")
s2.ABC_rej <- ABC_rejection(
    model=s2.model,
    prior=s2.prior,
    nb_simul=nb_simul,
    summary_stat_target=s2.sum_stat_obs,
    tol=tr)

print("Fine simulazione ABC 2")
s2.ABC_param = data.frame("p1" = s2.ABC_rej$param[,1], "p2" = s2.ABC_rej$param[,2])

p=ggplot(s2.ABC_param, aes(x=p1))+
    
    geom_histogram(color="darkblue", fill="lightblue") +
    
    geom_vline(aes(xintercept = mean(p1)), color = "red", linetype="longdash") +
    
    geom_vline(aes(xintercept = mean(s2.gr)), color = "blue", linetype="longdash") +
    
    labs(y = "counts", x = expression(theta[1])) +
    
    theme_minimal()

ggsave("estimated_p1_sim2.png", plot=p, path = paste0(path, "/plots"))

p=ggplot(s2.ABC_param, aes(x=p2)) +
    
    geom_histogram(color="darkblue", fill="lightblue") +
    
    geom_vline(aes(xintercept = mean(p2)), color = "red", linetype="longdash") +
    
    geom_vline(aes(xintercept = mean(s2.dr)), color = "blue", linetype="longdash") +
    
    labs(y = "counts", x = expression(theta[2])) +
    
    theme_minimal()

ggsave("estimated_p2_sim2.png", plot=p, path = paste0(path, "/plots"))


s2.tab <- cbind(s2.ABC_rej$param[,1], s2.ABC_rej$param[,2])
s2.tab <- as.table(s2.tab)
write.table(s2.tab, file = "sim2_result.csv", sep = ",", col.names = NA,
            qmethod = "double")
print("Fine simulazione 2")

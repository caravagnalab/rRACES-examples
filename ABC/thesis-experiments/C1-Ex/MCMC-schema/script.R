library(EasyABC)
library(rRACES)
library(ggplot2)

path = getwd()
print(path)

set.seed(1)

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

# SIMULAZIONE 1

# Parameter of the model
s1.gr = 0.08567
s1.dr = 0
s1.t=150
# Parameter simulation
nb_simul=3000
n=10

s1.prior = list(c("unif",0,0.1))

###########################################################
# Create simulation
sim <- new(Simulation, 'sim')
sim$history_delta = 1

# Add the mutant A
sim$add_mutant('A', s1.gr, s1.dr)

# Place cell on the tissue
dim <- sim$get_tissue_size()
sim$place_cell('A', dim[1]/2, dim[2]/2)

# Run the simulation for other t time units
sim$run_up_to_time(s1.t)

# Summary statistic of the simulation (number of cells)
s1.sum_stat_obs = c(sim$get_counts()$counts[1])

###################################################
p = plot_tissue(sim)
ggsave("tissue_plot_sim1.png", plot=p, path = paste0(path, "/plots"))

P = plot_muller(sim)
ggsave("muller_plot_sim1.png", plot=p, path = paste0(path, "/plots"))

######################################################
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

# Create the gr_model
s1.model = create_model(s1.t)

print("Inizio simulazione ABC 1")
s1.ABC_mcmc <-ABC_mcmc(method="Marjoram", 
    model=s1.model, 
    prior=s1.prior, 
    summary_stat_target=s1.sum_stat_obs, 
    n_rec=n)

print("Fine simulazione ABC 1")

s1.ABC_param = data.frame("p1" = s1.ABC_mcmc$param)

p = ggplot(s1.ABC_param, aes(x=p1)) +
    
    geom_histogram(color="darkblue", fill="lightblue") +
    
    geom_vline(aes(xintercept = mean(p1)), color = "red", linetype="longdash") +
    
    geom_vline(aes(xintercept = mean(s1.gr)), color = "blue", linetype="longdash") +
    
    labs(y = "counts", x = expression(theta)) +
    
    theme_minimal()

ggsave("estimated_parameter_sim1.png", plot=p, path = paste0(path, "/plots"))


s1.tab <- cbind(s1.ABC_mcmc$param[,1])
s1.tab <- as.table(s1.tab)
write.table(s1.tab, file = "sim1_result.csv", sep = ",", col.names = NA,
            qmethod = "double")


print("Fine simulazione 1")
#######################################################################################################
##############################################
#######################################################################################################
##############################################
#######################################################################################################

# SIMULAZIONE 2
set.seed(1)

print("Inizio simulazione 2")
# Parameter of the model
s2.gr = 0.0879
s2.dr = 0.00100546
s2.t=150
nb_simul=3000
n=10

###########################################################
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

# Summary statistic of the simulation (number of cells)
s2.sum_stat_obs = c(sim$get_counts()$counts[1])

###################################################
p = plot_tissue(sim)
ggsave("tissue_plot_sim2.png", plot=p, path = paste0(path, "/plots"))

p = plot_muller(sim)
ggsave("muller_plot_sim2.png", plot=p, path = paste0(path, "/plots"))

######################################################
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

# Create the gr_model
s2.model = create_model(s2.t)

# Let's assume that the parameter is in a uniform distribution [0,0.1]
s2.prior = list(c("unif",0,0.1), c("unif",0,0.05))

print("Inizio simulazione ABC 2")
s2.ABC_mcmc <- ABC_mcmc(method="Marjoram", 
    model=s2.model, 
    prior=s2.prior, 
    summary_stat_target=s2.sum_stat_obs, 
    n_rec=n)

print("Fine simulazione ABC 2")
s2.ABC_param = data.frame("p1" = s2.ABC_mcmc$param[,1], "p2" = s2.ABC_mcmc$param[,2])

p = ggplot(s2.ABC_param, aes(x=p1))+
    
    geom_histogram(color="darkblue", fill="lightblue") +
    
    geom_vline(aes(xintercept = mean(p1)), color = "red", linetype="longdash") +
    
    geom_vline(aes(xintercept = mean(s2.gr)), color = "blue", linetype="longdash") +
    
    labs(y = "counts", x = expression(theta[1])) +
    
    theme_minimal()

ggsave("estimated_p1_sim2.png", plot=p, path = paste0(path, "/plots"))

p = ggplot(s2.ABC_param, aes(x=p2)) +
    
    geom_histogram(color="darkblue", fill="lightblue") +
    
    geom_vline(aes(xintercept = mean(p2)), color = "red", linetype="longdash") +
    
    geom_vline(aes(xintercept = mean(s2.dr)), color = "blue", linetype="longdash") +
    
    labs(y = "counts", x = expression(theta[2])) +
    
    theme_minimal()

ggsave("estimated_p2_sim2.png", plot=p, path = paste0(path, "/plots"))


s2.tab <- cbind(s2.ABC_mcmc$param[,1], s2.ABC_mcmc$param[,2]) 
s2.tab <- as.table(s2.tab)
write.table(s2.tab, file = "sim2_result.csv", sep = ",", col.names = NA,
            qmethod = "double")

print("Fine simulazione 2")

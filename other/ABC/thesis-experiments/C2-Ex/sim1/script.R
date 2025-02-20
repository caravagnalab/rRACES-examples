library(EasyABC)
library(rRACES)
library(ggplot2)

# Setting working space -----------------------------------------------
path = getwd()
setwd(path)
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

print("Inizio simulazione")
# Simulations settings
set.seed(1)
gr1 = 0.021
gr2 = 0.0564
t1=150
t2=150

nb_simul=2000
tr=0.1
prior = list(c("unif", 50, 250),
             c("unif", 0.01, 0.07))


# Create simulation
sim <- new(Simulation, 'sim')
sim$history_delta = 1

# Add the mutant A
sim$add_mutant('A', gr1, 0)

# Place cell on the tissue
dim <- sim$get_tissue_size()
sim$place_cell('A', dim[1]/2, dim[2]/2)

# Run the simulation for other t time units
sim$run_up_to_time(sim$get_clock() + t1)

sim$add_mutant("B", gr2, 0)

sim$mutate_progeny(sim$choose_cell_in("A"), "B")

sim$run_up_to_time(sim$get_clock() + t2)
# Summary statistic of the simulation
# number of cells
sum_stat_obs = c(sim$get_counts()$counts[1],
                 sim$get_counts()$counts[2])
sum_stat_obs

p <- plot_tissue(sim)
ggsave("tissue_plot.png", plot = p, path = paste0(path, "/plots"))
print("Tissue plot saved")

p <- plot_muller(sim)
ggsave("muller_plot.png", plot = p, path = paste0(path, "/plots"))
print("Muller plot saved")

# Model definition ---------------------------------------------------
create_model <- function(gr, tf){
    model=function(x){
        if (all(x > 0)){
            tmp <- new(Simulation, 'tmp')
            tmp$add_mutant('A', gr, 0)
            
            dim <- tmp$get_tissue_size()
            tmp$place_cell("A", dim[1]/2, dim[2]/2)
            
            tmp$run_up_to_time(tmp$get_clock() + x[1])
            
            tmp$add_mutant("B", x[2], 0)
            
            tmp$mutate_progeny(tmp$choose_cell_in("A"), "B")
            
            tmp$run_up_to_time(tmp$get_clock() + (tf - x[1]))
            
            s=c(tmp$get_counts()$counts[1],
                tmp$get_counts()$counts[2])
        } else {
            s=c(0,0)
        }
    }
    print("The model has been created")
    return(model)
}


# Create the gr_model
m <- create_model(gr1, (t1+t2))

print("Inizio ABC rej")
ABC_rej <- ABC_rejection(
    model=m,
    prior=prior,
    nb_simul=nb_simul,
    summary_stat_target=sum_stat_obs,
    tol=tr)

print("Fine ABC rej")
# Plot the results
# Define the dataframe with the reuslt
param_df = data.frame("p1" = ABC_rej$param[,1], "p2" = ABC_rej$param[,2])

# P1 
p <- ggplot(param_df, aes(x = p1)) + 
    
    geom_histogram(color="darkblue", fill="lightblue") +
    
    geom_vline(aes(xintercept = t1), color = "blue", linetype = "longdash") + 
    
    geom_vline(aes(xintercept = mean(p1)), color = "red") +
    
    labs(y = "counts", x = expression(theta[1])) +
    
    #ggtitle(paste0("Accepted values for", expression(theta[1]))) +
    
    theme_minimal()

ggsave("p1_value_accepted.png", plot = p, path = paste0(path, "/plots"))

# P2
p <- ggplot(param_df, aes(x=p2)) + 
    
    geom_histogram(color="darkblue", fill="lightblue") +
    
    geom_vline(aes(xintercept = gr2), color = "blue", linetype = "longdash") + 
    
    geom_vline(aes(xintercept = mean(p2)), color = "red") +
    
    labs(y = "counts", x = expression(theta[2])) +
    
    #ggtitle(paste0("Accepted values for ", expression(theta[2]))) +
    
    theme_minimal()

ggsave("p2_value_accepted.png", plot = p, path = paste0(path, "/plots"))

tab <- cbind(param_df[,1], param_df[,2])
tab <- as.table(tab)
write.table(tab, file = "accepted_value.csv", sep = ",", col.names = NA,
            qmethod = "double")

print("Fine simulazione")

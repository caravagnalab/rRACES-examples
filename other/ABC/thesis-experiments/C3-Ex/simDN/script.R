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

# Paramaters settings
gr.A = 0.08
gr.B = 0.06

# first sample on the tissue (and t_therapy)
t1=150
# second sample on the tissue
t2=25
# add mutant B
t3=40
# third sample on the tissue
t4=50
# forth sample on the tissue
t5=150

# Create simulation
sim <- new(Simulation, 'sim')
sim$history_delta = 1

# Add the mutant A
sim$add_mutant('A', gr.A, 0)

# Place cell on the tissue
dim <- sim$get_tissue_size()
sim$place_cell('A', dim[1]/2, dim[2]/2)

# Evolve the tissue with mutant A until t1
sim$run_up_to_time(sim$get_clock() + t1)
p <- plot_tissue(sim)
ggsave("tissue_at_t1.png", plot = p, path = paste0(path, "/plots"))
print("saved tissue_at_t1.png")
# Register the first sample
sample.t1 <- sim$get_counts()$counts[1]

# start therapy
sim$update_rates("A", c(death=0.1))

# Evolve the tissue with mutant A for t2 time units
sim$run_up_to_time(sim$get_clock() + t2)
p <- plot_tissue(sim)
ggsave("tissue_at_t2.png", plot = p, path = paste0(path, "/plots"))
print("saved tissue_at_t2.png")

sample.t2 <- sim$get_counts()$counts[1]
sum_stat_1 <- (sample.t2-sample.t1)/sample.t1

# Evolve the tissue with mutant A for t2 time units
sim$run_up_to_time(sim$get_clock() + t3)
p <- plot_tissue(sim)
ggsave("tissue_at_t3.png", plot = p, path = paste0(path, "/plots"))
print("saved tissue_at_t3.png")

sim$add_mutant("B", gr.B, 0)
sim$mutate_progeny(sim$choose_cell_in("A"), "B")

sim$run_up_to_time(sim$get_clock() + t4)
p <- plot_tissue(sim)
ggsave("tissue_at_t4.png", plot = p, path = paste0(path, "/plots"))
print("saved tissue_at_t4.png")

sample.t4.A <- sim$get_counts()$counts[1]
sample.t4.B <- sim$get_counts()$counts[2]

sum_stat_2 <- (sample.t4.A-sum_stat_1)/sample.t4.A
sum_stat_3 <- sample.t4.B

sim$run_up_to_time(sim$get_clock() + t5)
p <- plot_tissue(sim)
ggsave("tissue_at_t5.png", plot = p, path = paste0(path, "/plots"))
print("saved tissue_at_t5.png")

sample.t5 <- sim$get_counts()$counts[2]
sum_stat_4 <- sample.t5

p <- plot_muller(sim)
ggsave("muller.png", plot = p, path = paste0(path, "/plots"))
print("saved muller.png")

sum_stat_obs = c(sum_stat_1, sum_stat_2, sum_stat_3, sum_stat_4)

# Model definition ---------------------------------------------------
create_model <- function(x1, y1, y2, y4, y5){
    model=function(x){
        if (all(x > 0)){
            # Create simulation
            tmp <- new(Simulation, 'tmp')
            # Add the mutant A
            tmp$add_mutant('A', x1, 0)
            
            # Place cell on the tissue
            dim <- tmp$get_tissue_size()
            tmp$place_cell('A', dim[1]/2, dim[2]/2)
            
            # Evolve the tissue with mutant A until t1
            tmp$run_up_to_time(tmp$get_clock() + y1)
            # Register the first sample
            sample.t1 <- tmp$get_counts()$counts[1]
            
            # start therapy
            tmp$update_rates("A", c(death=0.1))
            
            # Evolve the tissue with mutant A for t2 time units
            tmp$run_up_to_time(tmp$get_clock() + y2)
            
            sample.t2 <- tmp$get_counts()$counts[1]
            sum_stat_1 <- (sample.t2-sample.t1)/sample.t1
            
            # Evolve the tissue with mutant A for t3 time units
            tmp$run_up_to_time(tmp$get_clock() + (x[1]-(t1+t2)))
            
            tmp$add_mutant("B", x[2], 0)
            tmp$mutate_progeny(tmp$choose_cell_in("A"), "B")
            
            tmp$run_up_to_time(tmp$get_clock() + y4)
            
            sample.t4.A <- tmp$get_counts()$counts[1]
            sample.t4.B <- tmp$get_counts()$counts[2]
            
            sum_stat_2 <- (sample.t4.A-sum_stat_1)/sample.t4.A
            sum_stat_3 <- sample.t4.B
            
            tmp$run_up_to_time(tmp$get_clock() + y5)
            
            sample.t5 <- tmp$get_counts()$counts[2]
            sum_stat_4 <- sample.t5
            
            sum_stat_obs = c(sum_stat_1, sum_stat_2, sum_stat_3, sum_stat_4)
        } else {
            s=c(0,0,0,0)
        }
    }
    print("The model has been created")
    return(model)
}

# Create the gr_model
m <- create_model(gr.A, t1, t2, t4, t5)

nb_simul=2000
tr=0.1
prior = list(c("unif", 150, 400), c("unif", 0.03, 0.1))

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
    
    geom_vline(aes(xintercept = (t1+t2+t3)), color = "blue", linetype = "longdash") +
    
    geom_vline(aes(xintercept = mean(p1)), color = "red") +
    
    labs(y = "counts", x = expression(theta[1])) +
    
    #ggtitle(paste0("Accepted values for", expression(theta[1]))) +
    
    theme_minimal()

ggsave("t3_value_accepted.png", plot = p, path = paste0(path, "/plots"))

# P2
p <- ggplot(param_df, aes(x = p2)) +
    
    geom_histogram(color="darkblue", fill="lightblue") +
    
    geom_vline(aes(xintercept = gr.B), color = "blue", linetype = "longdash") +
    
    geom_vline(aes(xintercept = mean(p2)), color = "red") +
    
    labs(y = "counts", x = expression(theta[2])) +
    
    #ggtitle(paste0("Accepted values for", expression(theta[1]))) +
    
    theme_minimal()

ggsave("grB_value_accepted.png", plot = p, path = paste0(path, "/plots"))

tab <- cbind(param_df[,1], param_df[,2])
tab <- as.table(tab)
write.table(tab, file = "accepted_value.csv", sep = ",", col.names = NA,
            qmethod = "double")

print("Fine simulazione")

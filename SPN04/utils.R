
run_up_to_size_by_steps <- function(sim, mutant_name, size, delta_time = 1) {
  while ((sim$get_cells() %>% dplyr::filter(mutant == mutant_name) %>% nrow()) < size) {
    sim$run_up_to_time(sim$get_clock() + delta_time)
  }
  sim
}

run_down_to_size_by_steps <- function(sim, mutant_name, size, delta_time = 1) {
  while ((sim$get_cells() %>% dplyr::filter(mutant == mutant_name) %>% nrow()) > size) {
    sim$run_up_to_time(sim$get_clock() + delta_time)
  }
  sim
}

# Sweep function ####
sweep_population <- function(sim, old_mutant, new_mutant, desired_fraction_sweep, first_reduction = 2, delta_time = 1) {
  old_gr <- sim$get_species() %>% dplyr::filter(mutant == old_mutant) %>% pull(growth_rate)
  sim$update_rates(old_mutant, rates = c(growth = old_gr / first_reduction))
  new_counts <- sim$get_counts() %>% dplyr::filter(mutant == new_mutant) %>% pull(counts)
  total_counts <- sim$get_counts() %>% pull(counts) %>% sum()
  frac_new <- new_counts / total_counts
  while (frac_new <= desired_fraction_sweep) {
    sim$run_up_to_time(sim$get_clock() + delta_time)
    
    old_gr <- sim$get_species() %>% dplyr::filter(mutant == old_mutant) %>% pull(growth_rate)
    new_counts <- sim$get_counts() %>% dplyr::filter(mutant == new_mutant) %>% pull(counts)
    total_counts <- sim$get_counts() %>% pull(counts) %>% sum()
    frac_new <- new_counts / total_counts
    sim$update_rates(old_mutant, rates = c(growth = old_gr * (1 - frac_new)))
  }
  sim$update_rates(old_mutant, rates = c(growth = 0, death=5))
  sim
}
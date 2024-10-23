#This script is for simulating a lot of infections
#based on different connectance, spatial coverage, 
#and impact on primary.


###THE PARAMETERS THAT I'M INTERESTED IN
impact_on_primary <- seq(0,1,0.10)
spatial_coverage <- seq(0,1,0.10)

expanded_param <- expand.grid(
  impact_P = impact_on_primary, 
  sp_cov = spatial_coverage,
  connectance = seq(1,10),
  replicate = seq(1,10))

###The function to simulate a lot of infection
simulate_infection <- function(params) {
  impact_P <- params$impact_P
  sp_cov <- params$sp_cov
  connectance <- params$connectance
  replicate <- params$replicate
  
  # Call your discrete_trito_model with the correct parameters
  return(
    discrete_trito_model(1000, param, 500, impact_P, 0.9, sp_cov, 
                              igraph_list[[connectance]][[replicate]], 1)
    )
}

#Batch processing
for (batch in 1:num_batches) {
  # Determine start and end indices for the current batch
  start <- (batch - 1) * batch_size + 1
  end <- min(batch * batch_size, nrow(expanded_param))
  
  # Run simulations in parallel for the current batch
  results_batch <- mclapply(start:end, function(i) {
    simulate_infection(expanded_param[i, ])
  }, mc.cores = num_cores - 2)
  
  # Save the results of the current batch
  save(results_batch, 
       file = here("Output", "Connectivity_Intensity", 
                   paste0("sim_result_", batch, ".RData")))
  
  # Remove the current batch results from memory
  rm(results_batch)
}


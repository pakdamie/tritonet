# --------------------------------------------------------------------------------
# Research Question:
# How does the disturbance intensity of the primary vector influence the effective
# reproductive number (RE)? Additionally, how does this change based on the
# transmission efficiency of the secondary vector?
# --------------------------------------------------------------------------------

# Define the different parameter scenarios for simulation
param_interest <- c("standard", "nonesec", "worse_m", "no_diff", "better_m")

# Generate a list of parameter sets based on the specified scenarios
param_list <- lapply(param_interest, function(x) get_parameters(x))

# Define different mortality rates for the primary vector
Mortality_P <- data.frame(Mort_P = c(0, 0.01, 0.25, 0.5, 0.75))

# When is the system disturbed?
dstb_time <- get_parameters("standard")[["disturbance_time"]]

# --------------------------------------------------------------------------------
# Step 1: Simulate Model Output and Calculate Re
# --------------------------------------------------------------------------------

# Iterate over each parameter set and simulate model output.
# Then, compute the expanded human effective reproductive number (Re).
RE_mortality_P <- lapply(param_list, function(x) {
  # Simulate model output for a given parameter set and mortality rate
  sim_output <- Simulate_Model_Output(x,
    infection_start = "No", # No initial infection
    variable_interest = "mortality_P", # Focus on primary vector mortality
    vector_value = Mortality_P # Use defined mortality rates
  )

  # Compute the Human RE from the simulated output
  lapply(sim_output, Calculate_Human_Reff_Expanded, param = x)
})

# --------------------------------------------------------------------------------
# Step 2: Assign Parameter Type and Mortality Rate to Each List Element
# --------------------------------------------------------------------------------

# Loop through each parameter set
for (i in seq_along(param_interest)) {
  for (j in seq_len(nrow(Mortality_P))) {
    # Assign the parameter type (e.g., "standard", "nonesec", etc.)
    RE_mortality_P[[i]][[j]]$param <- param_interest[i]

    # Assign the corresponding mortality rate value
    RE_mortality_P[[i]][[j]]$id <- Mortality_P[j, ]
  }
}

# Convert the nested list structure into a data frame
RE_mortality_P_DF <- do.call(rbind, unlist(RE_mortality_P, recursive = FALSE))

# --------------------------------------------------------------------------------
# Step 3: Subset Time Period of Interest
# --------------------------------------------------------------------------------

# Extract a specific time window from the simulation results.
# The subset includes data slightly before and after the disturbance event.
RE_mortality_P_post <- subset(
  RE_mortality_P_DF,
  RE_mortality_P_DF$time > dstb_time - 100 & # Start time
    RE_mortality_P_DF$time < (dstb_time + 300) # End time (disturbance + 300 days)
)

# --------------------------------------------------------------------------------
# Step 4: Save Processed Data
# --------------------------------------------------------------------------------

# Save the filtered dataset as an RDS file in the "Output" directory

saveRDS(RE_mortality_P_post, file = here("Output", "df_RE_mortality_P_R0.rds"))

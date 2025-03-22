# --------------------------------------------------------------------------------
# Research Question:
# How does the secondary contribution to the RE change over time and
# with the change in the transmission efficiency?
# --------------------------------------------------------------------------------
param_standard <- get_parameters("standard")
Mortality_P <- seq(0.01, 0.99, length = 25)
modifier <- seq(0.05, 2, 0.01)

### Pull out the standard parameters
f_M_standard <- param_standard["f_M"] # Biting rate of s. vector
theta_M_standard <- param_standard["theta_M"] # Trans. prob of s.vector
f_P_standard <- param_standard["f_P"] # Biting rate of p. vector
theta_P_standard <- param_standard["theta_P"] # Trans. prob of p.vector

# A way of figuring out what parameters to vary
f_M_interest <- ((f_P_standard * theta_P_standard) * modifier) / theta_M_standard

secondary_param <-
  data.frame(expand.grid(
    f_M = f_M_interest,
    mortality_P = Mortality_P
  ))

secondary_param_list <- vary_parameter_value(
  param_standard,
  c("f_M", "mortality_P"), secondary_param
)

# --------------------------------------------------------------------------------
# Step 1: Simulate Model Output and Calculate Re
# --------------------------------------------------------------------------------

RE_SM <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("f_M", "mortality_P"),
    vector_value = secondary_param
  ) |>
  Calculate_change_baseline(
    secondary_param_list,
    secondary_param, "No"
  )
saveRDS(RE_SM, file = here("Output", "RE_SM.rds"))


# --------------------------------------------------------------------------------
# Step 2: Simulate Model Output (Temporal calculation of RE)
# --------------------------------------------------------------------------------


RE_SM_2 <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("f_M", "mortality_P"),
    vector_value = secondary_param
  ) |>
  Calculate_change_baseline(
    secondary_param_list,
    secondary_param, "Yes"
  )

# We are getting the standardized ratio that we want
RE_SM_2$standardized_ratio <- (RE_SM_2$f_M * theta_M_standard) / (f_P_standard * theta_P_standard)
RE_SM_2_subset <- subset(
  RE_SM_2,
  RE_SM_2$mortality_P %in% c(0.01) &
    RE_SM_2$standardized_ratio %in%
      c(0.05, 0.5, 1, 1.5, 2)
)

saveRDS(RE_SM_2_subset, file = here("Output", "RE_SM_2_subset.rds"))

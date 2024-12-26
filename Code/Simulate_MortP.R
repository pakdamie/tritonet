# How does the disturbance intensity of the primary vector influence RE?

# Retrieve "standard" parameters and set the different values of disturbances
param_standard <- get_parameters("standard")
Mortality_P <- data.frame(c(0.01, seq(0.05, 1, 0.05)))

# Simulate model output and calculate the RE
RE_mortality_P <-
  Simulate_Model_Output(param_standard, "mortality_P", Mortality_P) |>
  lapply(Calculate_Human_REff, param = param_standard)

# Assign mortality_P value to each list element
for (i in seq_along(RE_mortality_P)) {
  RE_mortality_P[[i]]$id <- Mortality_P[i, ]
}
RE_mortality_P_DF <- do.call(rbind, RE_mortality_P)

### Subset time period of interest (slightly before disturbance and post disturbance)
RE_mortality_P_post <- subset(
  RE_mortality_P_DF,
  RE_mortality_P_DF$time > 9000 &
    RE_mortality_P_DF$time < 14000
)

RE_mortality_P_post$NewHumanCases <-
  ((param_standard["f_P"] * param_standard["theta_P"] *
    (RE_mortality_P_post$NP - RE_mortality_P_post$PS) / 1000) +
    (param_standard["f_M"] * param_standard["theta_M"] *
      (RE_mortality_P_post$NM - RE_mortality_P_post$MS) / 1000)) * RE_mortality_P_post$HS

saveRDS(RE_mortality_P_post, file = here("Output", "RE_mortality_P_post.rds"))


## What if the parameters are the same for both the primary and secondary vectors?

param_same <- get_parameters("no_diff")

RE_mortality_P_same <-
  Simulate_Model_Output_PostD( "mortality_P", Mortality_P) |>
  lapply(Calculate_Human_REff, param =  get_parameters("post_disturb"))


# Assign mortality_P value to each list element
for (i in seq_along(RE_mortality_P_same)) {
  RE_mortality_P_same[[i]]$id <- Mortality_P[i, ]
}
RE_mortality_P_same_DF <- do.call(rbind, RE_mortality_P_same)

RE_mortality_P_same_DF <- subset(
        RE_mortality_P_same_DF,
        RE_mortality_P_same_DF$time >100 &
        RE_mortality_P_same_DF$time < 3000
)
RE_mortality_P_same_DF$NewHumanCases <-
  ((param_same["f_P"] * param_same["theta_P"] *
    (RE_mortality_P_same_DF$NP - RE_mortality_P_same_DF$PS) / 1000) +
    (param_same["f_M"] * param_same["theta_M"] *
      (RE_mortality_P_same_DF$NM - RE_mortality_P_same_DF$MS) / 1000)) * RE_mortality_P_same_DF$HS

saveRDS(RE_mortality_P_same_DF, file = "Output/RE_mortality_P_same_DF.rds")

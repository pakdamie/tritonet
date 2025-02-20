# How does the disturbance intensity of the primary vector influence RE?

# Retrieve "standard" parameters and set the different values of disturbances
param_standard <- get_parameters("standard")
Mortality_P <- data.frame(Mort_P = c(0, 0.01, 0.25, 0.5, 0.75))

# Simulate model output and calculate the RE from the output
RE_mortality_P <-
  Simulate_Model_Output(param_standard,
    infection_start = "No",
    "mortality_P",
    Mortality_P) |>
  lapply(Calculate_Human_Reff_Expanded, param = param_standard)

# Assign mortality_P value to each list element
for (i in seq_along(RE_mortality_P)) {
  RE_mortality_P[[i]]$id <- Mortality_P[i, ]
}
RE_mortality_P_DF <- do.call(rbind, RE_mortality_P)

### Subset time period of interest (slightly before disturbance and post disturbance)
RE_mortality_P_post <- subset(
  RE_mortality_P_DF,
  RE_mortality_P_DF$time > 9100 &
  RE_mortality_P_DF$time < (9125 + 300)
)


saveRDS(file = here("Output", "df_RE_mortality_P_R0.rds"), RE_mortality_P_post)

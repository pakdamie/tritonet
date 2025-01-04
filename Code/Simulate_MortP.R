# How does the disturbance intensity of the primary vector influence RE?
# Standard case
# Retrieve "standard" parameters and set the different values of disturbances
param_standard <- get_parameters("standard")
Mortality_P <- data.frame(Mort_P = c(0, 0.10,  0.25, 0.50, 0.75))

#Simulate model output and calculate the RE
RE_mortality_P <-
  Simulate_Model_Output(param_standard, "mortality_P", Mortality_P) |>
  lapply(Calculate_Human_REff, param = param_standard) 
  
#Assign mortality_P value to each list element
for (i in seq_along(RE_mortality_P)) {
  RE_mortality_P[[i]]$id <- Mortality_P[i, ]
}
RE_mortality_P_DF <- do.call(rbind, RE_mortality_P)

### Subset time period of interest (slightly before disturbance and post disturbance)
RE_mortality_P_post <- subset(
  RE_mortality_P_DF,
  RE_mortality_P_DF$time > 8900 &
  RE_mortality_P_DF$time <  9125 + 2000
)

saveRDS(RE_mortality_P_post, file = here("Output", "RE_mortality_P_post.rds"))


## What if the parameters are the same for both the primary and secondary vectors?

RE_mortality_P_same <-
  Simulate_Model_Output_PostD("mortality_P", Mortality_P) |>
  lapply(Calculate_Human_REff, param = get_parameters("post_disturb"))


# Assign mortality_P value to each list element
for (i in seq_along(RE_mortality_P_same)) {
  RE_mortality_P_same[[i]]$id <- Mortality_P[i, ]
}
RE_mortality_P_same_DF <- do.call(rbind, RE_mortality_P_same)

RE_mortality_P_same_DF <- subset(
 RE_mortality_P_same_DF,
 RE_mortality_P_same_DF$time > 100 &
 RE_mortality_P_same_DF$time < 3000
)

## What if the secondary vectors do not contribute at all to disease?

Mortality_ThetaM = Mortality_P
Mortality_ThetaM$theta_M = 0

RE_mortality_P_same <-
  Simulate_Model_Output_PostD(c("mortality_P", "theta_M"), Mortality_ThetaM) |>
  lapply(Calculate_Human_REff, param = get_parameters("post_disturb"))




saveRDS(RE_mortality_P_same_DF, file = "Output/RE_mortality_P_same_DF.rds")

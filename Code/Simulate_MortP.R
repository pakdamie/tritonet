# How does the mortality of the primary vector influence RE?

#Get the standard parameters
param_standard <- get_parameters("standard")
Mortality_P <- data.frame(c(0.01, seq(0.05, 1, 0.05)))

###Simulate model output
Mod_Mort_P <- Simulate_Model_Output(param_standard, "mortality_P", Mortality_P)

###Calculate human RE
RE_mortality_P <- lapply(Mod_Mort_P, 
                         Calculate_Human_REff, 
                         param = param_standard)

# Assign mortality_P value
for (i in seq_along(RE_mortality_P)) {
  RE_mortality_P[[i]]$id <- Mortality_P[i,]
}

RE_mortality_P_DF <- do.call(rbind, RE_mortality_P)

###Subset time period of interest
RE_mortality_P_post <- subset(
  RE_mortality_P_DF,
  RE_mortality_P_DF$time > 9000 &
  RE_mortality_P_DF$time < 14000)

RE_mortality_P_post$NewHumanCases <- 
  ((param_standard["f_P"] * param_standard["theta_P"] * 
   (RE_mortality_P_post$NP - RE_mortality_P_post$PS)/1000) +   
  (param_standard["f_M"] * param_standard["theta_M"] * 
   (RE_mortality_P_post$NM - RE_mortality_P_post$MS)/1000)) * RE_mortality_P_post$HS

saveRDS(RE_mortality_P_post, file = here("Output","RE_mortality_P_post.rds"))

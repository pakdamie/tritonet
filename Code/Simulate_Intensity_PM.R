### Simulate varying intensity of disturbances for both primary and secondary

param_standard <- get_parameters("standard")
Mortality_P <- c(seq(0, 1, 0.1))
Mortality_M <- c(seq(0, 1, 0.1))
Mortality_PM_Grid <- expand.grid(Mortality_P, Mortality_M)

### Simulate model output
Mod_Mort_PM <- Simulate_Model_Output(
  param_standard,infection_start = "No", c("mortality_P", "mortality_M"), 
  Mortality_PM_Grid
)

### Calculate human RE
RE_Mort_PM <- lapply(Mod_Mort_PM, Calculate_Human_Reff_Expanded, param = param_standard)

# Assign mortality_P value
for (i in 1:nrow(Mortality_PM_Grid)) {
  RE_Mort_PM[[i]]$Mort_P <- Mortality_PM_Grid[i, 1]
  RE_Mort_PM[[i]]$Mort_M <- Mortality_PM_Grid[i, 2]
}
RE_Mort_PM_DF <- do.call(rbind, RE_Mort_PM)

### Subset time period of interest
RE_Mort_PM_DF <- subset(
  RE_Mort_PM_DF,
  RE_Mort_PM_DF$time > 9100 &
    RE_Mort_PM_DF$time < 14000
)

Maximum_RE_Mort_PM <- do.call(
  rbind,
  by(RE_Mort_PM_DF,
    RE_Mort_PM_DF[, c("Mort_P", "Mort_M")],
    function(x) {
      x[which.max(x$RE), ]
    },
    simplify = TRUE
  )
)


ggplot(Maximum_RE_Mort_PM, aes(x = (1-Mort_P), y= (1-Mort_M), fill = RE)) + 
  geom_tile() + 
  scale_fill_viridis()



save(Maximum_RE_Mort_PM, file = "Output/Maximum_RE_Mort_PM.rds")

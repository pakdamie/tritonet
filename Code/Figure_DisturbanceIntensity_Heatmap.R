###Figure_2_Varying intensity of disturbances

# Running models and cleaning it up ---------------------------------------
param_standard<- get_parameters()
Mortality_P <- c(0.01, seq(0.05, 1, 0.05))
Mortality_M <- c(0.01, seq(0.05, 1, 0.05))
Mortality_PM_Grid <- expand.grid(Mortality_P, Mortality_M)

###Simulate model output
Mod_Mort_PM <- Simulate_Model_Output(
  param_standard, c("mortality_P", "mortality_M"), Mortality_PM_Grid)

###Calculate human RE
RE_Mort_PM <- lapply(Mod_Mort_PM, Calculate_Human_REff, param = param_standard)

# Assign mortality_P value
for (i in 1:nrow(Mortality_PM_Grid)) {
  RE_Mort_PM[[i]]$Mort_P <- Mortality_PM_Grid[i,1]
  RE_Mort_PM[[i]]$Mort_M <- Mortality_PM_Grid[i,2]
  
}
RE_Mort_PM_DF <- do.call(rbind, RE_Mort_PM)

###Subset time period of interest
RE_Mort_PM_DF  <- subset(
        RE_Mort_PM_DF,
        RE_Mort_PM_DF$time > 9120 &
        RE_Mort_PM_DF$time < 14000)

find_max_RE(RE_Mort_PM_DF )

Maximum_RE_Mort_PM <- do.call(
  rbind,
  by(RE_Mort_PM_DF,
    RE_Mort_PM_DF[,c("Mort_P", "Mort_M")],
    function(x) {
      x[which.max(x$RE), ]
    },
    simplify = TRUE
  )
)
   

# Figure ------------------------------------------------------------------

GG_Disturbance_Heatmap_RE <- 
  ggplot(Maximum_RE_Mort_PM,  
   aes(x = as.factor(Mort_P), 
       y = as.factor(Mort_M), 
       fill = RE)) + 
  geom_raster() + 
  scale_fill_viridis(option = 'rocket', name = expression(R[E])) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + 
  xlab(expression("Disturbance intensity of primary vector " * "(" *mu[P]* ")")) + 
  ylab(expression("Disturbance intensity of secondary vector " * "(" *mu[M]* ")")) + 
  coord_equal() + 
  theme(axis.text = element_text(color = 'black',size = 13),
        axis.title = element_text(color = 'black', size = 14));GG_Disturbance_Heatmap_RE

ggsave(here("Figures_Process", "Dec_10", "GG_Disturbance_Heatmap_RE.pdf"),
       width = 10, height = 10, units = 'in')


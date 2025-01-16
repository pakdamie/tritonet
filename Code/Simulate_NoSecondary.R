# How does the disturbance intensity of the primary vector influence RE?
# The secondary vector does not transmit
# Retrieve "standard" parameters and set the different values of disturbances
param_nonesec <- get_parameters("nonesec")
Mortality_P <- data.frame(Mort_P = c(0,0.01, 0.25, 0.5, 0.75))

#Simulate model output and calculate the RE
RE_mortality_P_ALT <-
        Simulate_Model_Output(param_nonesec , "mortality_P", Mortality_P) |>
        lapply(Calculate_Human_REff, param = param_nonesec ) 

#Assign mortality_P value to each list element
for (i in seq_along(RE_mortality_P_ALT)) {
        RE_mortality_P_ALT[[i]]$id <- Mortality_P[i, ]
}
RE_mortality_P_ALT_DF <- do.call(rbind, RE_mortality_P_ALT)

### Subset time period of interest (slightly before disturbance and post disturbance)
RE_mortality_P_ALT_post <- subset(
        RE_mortality_P_ALT_DF,
        RE_mortality_P_ALT_DF$time > 8900 &
                RE_mortality_P_ALT_DF$time <  9125 + 2000
)


RE_limits <- c(round(min(RE_mortality_P_ALT_post$RE),1), 
               round(max(RE_mortality_P_ALT_post$RE),1) + 0.1)

#Plot total vector abundance over time with the RE As color
Panel_A <- plot_NV_RE(RE_mortality_P_ALT_post,"No", RE_limits)

#Remove the situation where 100% of the primary vector is removed,
#not interesting (dynamically, and a little bug that occurs in Panel B! 
#If you do it individually, it works :/)
removed_0 <- subset(RE_mortality_P_ALT_post , RE_mortality_P_ALT_post $id != 0)

#Plot the secondary versus primary vector with the RE as color
Panel_B <- plot_NP_NM_RE(removed_0 , postdisturb = "No",RE_limits);

#Plot the secondary contribution to RE
Panel_C <- plot_NM_REff(removed_0, postdisturb = "No");

#Full composite figure with all
Panel_A + (Panel_B) + 
        plot_layout(guides='collect')  &
        theme(legend.position='right') 

ggsave("")
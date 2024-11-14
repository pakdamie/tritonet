# This is an analysis script for the one patch dynamics
###ONE PATCH SYSTEM

param_standard <-  
  c(b_H = 1/(1000), ##Human mortality rate
    b_P = 0.01, # P. Vector birth rate
    b_M = 0.01,  # S. Vector birth rate
    mu_H =1/(1000),  ##Human death rate
    mu_V = 0,  ##P. vector death rate
    f_P = 0.050, # Biting rate of the p. vector
    f_M = 0.040, # Biting rate of the s.vector
    theta_P = 0.90, # Transmission probability of p. vector
    theta_M = 0.90, # Transmission probability of s. vector
    theta_H  = 0.50, # Transmission probability of human
    gamma = 1/100,  # Recovery rate of infected human
    c_PM = 5e-6, ## Competition effect of p.vector on s.vector
    c_MP = 1e-6, ## Competition effect of s.vector on p.vector
    c_PP = 7e-6, ## Competition effect of p.vector on s.vector
    c_MM = 4e-6, ## Competition effect of s.vector on s.vector
    ntime = 10000,
    disturbance_time = 5000,
    delta_T = 1,
    prop = 1,
    mortality_P = 0.10,
    mortality_M = 0.95) 

#The survivorship of the 
mortality_P <- seq(0.1,1,0.10)

###Creating a list of parameters that I can loop through
parameter_mortality_P_list <- NULL

for(i in seq(1,length(mortality_P))){
   copied_param <- param_standard     
   copied_param['mortality_P']  <- mortality_P[i] 
   
   ## Put in list to be looped through   
   parameter_mortality_P_list[[i]] <- copied_param     
        
}

###Setting initial conditions
Initial_List <- create_initial_states(param_standard, patch_num = 1)

###Simulate the model with varying primary-removal 
model_output_mort_P <- NULL
for (p in seq(1,length(parameter_mortality_P_list))){

model_output_mort_P[[p]] <- 
  discrete_trito_model_rcpp_ONEPATCH(
        HS = Initial_List[[1]],
        HI = Initial_List[[2]],
        HR = Initial_List[[3]], 
        PS = Initial_List[[4]],
        PI = Initial_List[[5]],
        MS = Initial_List[[6]],
        MI = Initial_List[[7]],
        param = parameter_mortality_P_list[[p]])    

}

## Calculate the RE

RE_onepatch_List <- NULL

for (p in seq(1,length(parameter_mortality_P_list))){
        
  RE_onepatch_List [[p]] <- cbind.data.frame(
    calculate_R_effective_modified(
    model_output_mort_P[[p]],
    parameter_mortality_P_list[[p]]),
    primary_removal = mortality_P[[p]])
}

RE_onepatch_DF <- do.call(rbind,RE_onepatch_List )

# Figure showing how the RE changes over time with 
# changing total secondary vectors primary vectors
ggplot(subset(RE_onepatch_DF,
              RE_onepatch_DF$time > 4999),
   aes(x = N_P, y = N_M, color = RE, 
       group = primary_removal)) + 
   geom_path(size = 1.1) + 
   geom_point(data=
     subset(RE_onepatch_DF, 
            RE_onepatch_DF$time == 5000),
       aes(x = N_P, y = N_M), 
           size = 2.5, 
           color = 'black') + 
   scale_colour_gradient2(low = "#5a47eb", 
                          mid = 'lightgrey' ,
                          high = "#e72378", 
                          midpoint =1 ) + 
   theme_classic() + 
   xlab("Primary vectors") + 
   ylab("Secondary vectors") + 
   theme(legend.position = 'top') 


## At what abundance is RE the maximum?
eq_RE_onepatch <- subset(RE_onepatch_DF,
                         RE_onepatch_DF$time > 4999)

eq_RE_onepatch_df <- do.call(
  rbind,
  lapply(split(eq_RE_onepatch, 
    eq_RE_onepatch$primary_removal),
    function(x) x[which.max(x$RE),]))

ggplot(eq_RE_onepatch_df, 
    aes(x = N_P, y = N_M, fill = RE)) + 
  geom_point(size = 4, shape = 21, color = 'black') + 
  xlab("Total primary vectors") + 
  ylab("Total secondary vectors") + 
  scale_fill_viridis(option = 'rocket') + 
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15))
  
### Time
ggplot(eq_RE_onepatch_df, 
       aes(x = primary_removal, y = time - 5000, fill = RE)) + 
  geom_point(size = 4, shape = 21, color = 'black') + 
  xlab("Proportion of primary surviving") + 
  ylab("Time after disturbance") + 
  scale_fill_viridis(option = 'rocket') + 
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15)) + 
  ggtitle("When is the peak RE reached?")


# Standard Parameter ---------------------------------------------------------------
param_standard_2 <- c(
        b_H = 1 / (1000), ## Human mortality rate
        b_P = 0.01, # P. Vector birth rate
        b_M = 0.01, # S. Vector birth rate
        mu_H = 1 / (1000), ## Human death rate
        f_P = 0.02, # Biting rate of the p. vector
        f_M = 0.02 , # Biting rate of the s.vector
        theta_P = 0.70, # Transmission probability of p. vector
        theta_M = 0.70, # Transmission probability of s. vector
        theta_H = 0.50, # Transmission probability of human
        gamma = 1 / 90, # Recovery rate of infected human
        c_PM = 4e-6, ## Competition effect of p.vector on s.vector
        c_MP = 2e-6, ## Competition effect of s.vector on p.vector
        c_PP = 4.5e-6, ## Competition effect of p.vector on s.vector
        c_MM = 3e-6, ## Competition effect of s.vector on s.vector
        ntime = 365 * 50,
        disturbance_time = 365 * 25,
        delta_T = 1,
        prop = 1,
        mortality_P = 0.25, # This will change
        mortality_M = 1
)

# Running models and cleaning it up ---------------------------------------
modifier <- seq(0.10, 2, 0.20)
modifier_full <- as.numeric(modifier * param_standard_2[["theta_M"]])
Mortality_P <- seq(0, 1, 0.25)
Mortality_P[1] <- 0.01
theta_M_disturbance <- expand.grid(theta_M = modifier_full, 
                                   mortality_P = Mortality_P)


#Initial states
Initial_List <- create_initial_states(param_standard_2)






RE_theta_M_P_List <- NULL

for (param in seq(1:nrow(theta_M_disturbance))){
  
  param_changed <- param_standard_2   
  param_changed["theta_M"] <- theta_M_disturbance[param,]$theta_M
  param_changed["mortality_P"] <- theta_M_disturbance[param,]$mortality_P      
        
  model_f_M_disturbance <- 
    discrete_trito_model_rcpp_ONEPATCH(
     HS = Initial_List[[1]],
     HI = Initial_List[[2]],
     HR = Initial_List[[3]],
     PS = Initial_List[[4]],
     PI = Initial_List[[5]],
     MS = Initial_List[[6]],
     MI = Initial_List[[7]],
     param = param_changed)
  
  plot_list_groups(model_f_M_disturbance)
        
model_RE <- Calculate_analytical_REff(model_f_M_disturbance , param_changed)

model_RE$theta_M <- theta_M_disturbance[param,]$theta_M  
model_RE$mortality_P <- theta_M_disturbance[param,]$mortality_P  

RE_theta_M_P_List[[param]] = model_RE

}

RE_theta_M_P_DF <- do.call(rbind, RE_theta_M_P_List)


RE_theta_M_P_DF  <- subset(
        RE_theta_M_P_DF ,
        RE_theta_M_P_DF $time > 9124 &
        RE_theta_M_P_DF $time < 14000
)


RE_THETA_M_MAX_RE <- do.call(rbind,lapply(split(RE_theta_M_P_DF , 
                list(RE_theta_M_P_DF$theta_M, RE_theta_M_P_DF $mortality_P), drop = TRUE),
       function(x)
        x[which.max(x$RE),]))

ggplot(RE_THETA_M_MAX_RE, 
  aes(x = theta_M, y= as.factor(mortality_P), fill = RE^2)) + geom_tile() + 
  xlab("Secondary vector's competence compared to primary vector") + 
  ylab("Disturbance intensity of primary vector") + 
  scale_fill_gradient2(
          low = "#448C81",
          mid = "#EFEAD7",
          high = "#C25D31",
          midpoint = 1,
          name = expression(R[E])
  ) +
 scale_y_discrete(expand= c(0,0)) + 
 theme(axis.text = element_text(size = 14))
  

ggplot(RE_theta_M_P_DF ,aes(x = time, y= RE, group = as.factor(theta_M), 
                            color = as.factor(theta_M))) + 
        geom_path(size = 1.2) + facet_wrap(~mortality_P) +
        scale_color_viridis(discrete = TRUE)

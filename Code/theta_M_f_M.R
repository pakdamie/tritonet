
f_m_standard <- get_parameters("no_diff")["f_M"]
theta_m_standard <- get_parameters("no_diff")["theta_M"]


# Running models and cleaning it up ---------------------------------------
modifier_values <- seq(0.05, 2.5, 0.10)
modified_theta_M <- as.numeric(modifier_values * get_parameters("no_diff")[["theta_M"]])
modified_f_M <- as.numeric(modifier_values * get_parameters("no_diff")[["f_M"]])
Mortality_P <- 0.01
thetaf_M_disturbance <- expand.grid(theta_M = modified_theta_M, 
                                   f_M = modified_f_M,
                                   mortality_P = Mortality_P)

Mod_MortP_thetafM <- Simulate_Model_Output_PostD (
        c("theta_M", "f_M", "mortality_P"), 
        thetaf_M_disturbance)


RE_thetaf_M_P_List <- NULL

for (param in seq(1:nrow(thetaf_M_disturbance))){
        
  param_changed <- get_parameters("post_disturb")
  param_changed["theta_M"] <- thetaf_M_disturbance[param,]$theta_M
  param_changed["f_M"] <- thetaf_M_disturbance[param,]$f_M
  param_changed["mortality_P"] <- thetaf_M_disturbance[param,]$mortality_P   
  
 
  model_RE <- Calculate_Human_REff(Mod_MortP_thetafM[[param]], param_changed)
  
  model_RE$f_M <- thetaf_M_disturbance[param,]$f_M  
  model_RE$theta_M <- thetaf_M_disturbance[param,]$theta_M  
  model_RE$mortality_P <- thetaf_M_disturbance[param,]$mortality_P  
  
  RE_thetaf_M_P_List[[param]] = model_RE
        
}

RE_thetaf_M_P_DF <- do.call(rbind, RE_thetaf_M_P_List)







RE_THETAF_M_MAX_RE <- 
        do.call(rbind,lapply(split(RE_thetaf_M_P_DF , 
                                   list(RE_thetaf_M_P_DF$theta_M, 
                                        RE_thetaf_M_P_DF$f_M), 
                                   drop = TRUE),function(x)
                                           x[which.max(x$RE),]))
thetaF_M_contribution_to_RE_GG <- 
        ggplot(RE_THETAF_M_MAX_RE, 
      aes(x = as.factor(theta_M/theta_m_standard), 
          y= as.factor(f_M/f_m_standard ), 
          fill =RE)) + 
        geom_raster() + 
        xlab(expression("Multiplier of transmission probability " * "(" *theta[M]* ")")) + 
        ylab(expression("Multiplier of biting rate" * "(" *f[M]* ")")) + 
        scale_fill_viridis(option = 'rocket', name = "Secondary \ncontribution to RE") + 
        scale_y_discrete(expand= c(0,0)) + 
        scale_x_discrete(expand = c(0,0)) + 
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 15)); thetaF_M_contribution_to_RE_GG

x[which.max(x$RE),]))
eq_initial_list <- calculate_predisturb_initial_list()

# Running models and cleaning it up ---------------------------------------
modifier_values <- seq(0.05, 2.5, 0.10)
modified_theta_M <- as.numeric(modifier * get_parameters("no_diff")[["theta_M"]])
modified_f_M <- as.numeric(modifier * get_parameters("no_diff")[["f_M"]])
Mortality_P <- 0.10
thetaf_M_disturbance <- expand.grid(theta_M = modified_theta_M, 
                                   f_M = modified_f_M,
                                   mortality_P = Mortality_P)


Mod_MortP_thetafM <- Simulate_Model_Output(
        get_parameters("post_disturb"), 
        c("theta_M", "f_M", "mortality_P"), 
        thetaf_M_disturbance)


RE_thetaf_M_P_List <- NULL

for (param in seq(1:nrow(thetaf_M_disturbance))){
        
  param_changed <- get_parameters("post_disturb")
  param_changed["theta_M"] <- thetaf_M_disturbance[param,]$theta_M
  param_changed["f_M"] <- thetaf_M_disturbance[param,]$f_M
  param_changed["mortality_P"] <- thetaf_M_disturbance[param,]$mortality_P   
  
  Model_output<-
                discrete_trito_model_rcpp_ONEPATCH(
                        HS = empty_initial_list[[1]],
                        HI = empty_initial_list[[2]],
                        HR = empty_initial_list[[3]],
                        PS = empty_initial_list[[4]],
                        PI = empty_initial_list[[5]],
                        MS = empty_initial_list[[6]],
                        MI = empty_initial_list[[7]],
                        param =  param_changed)
        
        model_RE <- Calculate_Human_REff(Model_output, param_changed)
        
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
      aes(x = as.factor(theta_M/0.7), y= as.factor(f_M/0.02), 
          fill = MtoH/RE)) + 
        geom_raster() + 
        xlab(expression("Multiplier of transmission probability " * "(" *theta[M]* ")")) + 
        ylab(expression("Multiplier of biting rate" * "(" *f[M]* ")")) + 
        scale_fill_viridis(option = 'rocket', name = "Secondary \ncontribution to RE") + 
        scale_y_discrete(expand= c(0,0)) + 
        scale_x_discrete(expand = c(0,0)) + 
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 15)); thetaF_M_contribution_to_RE_GG

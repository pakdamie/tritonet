eq_initial_list <- calculate_predisturb_initial_list()

# Running models and cleaning it up ---------------------------------------
modifier_values <- seq(0.05, 2.5, 0.10)
modified_theta_M <- as.numeric(modifier * get_parameters("no_diff")[["theta_M"]])
Mortality_P <- seq(0, 1, 0.05)
Mortality_P[1] <- 0.01
theta_M_disturbance <- expand.grid(theta_M = modified_theta_M, 
                                   mortality_P = Mortality_P)


Mod_MortP_thetaM <- Simulate_Model_Output(
        get_parameters("post_disturb"), 
        c("theta_M", "mortality_P"), theta_M_disturbance)


RE_theta_M_P_List <- NULL

for (param in seq(1:nrow(theta_M_disturbance))){
        
        param_changed <- get_parameters("post_disturb")
        param_changed["theta_M"] <- theta_M_disturbance[param,]$theta_M
        param_changed["mortality_P"] <- theta_M_disturbance[param,]$mortality_P   
        
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
        
        model_RE$theta_M <- theta_M_disturbance[param,]$theta_M  
        model_RE$mortality_P <- theta_M_disturbance[param,]$mortality_P  
        
        RE_theta_M_P_List[[param]] = model_RE
        
}

RE_theta_M_P_DF <- do.call(rbind, RE_theta_M_P_List)


RE_THETA_M_MAX_RE <- 
        do.call(rbind,lapply(split(RE_theta_M_P_DF , 
                                   list(RE_theta_M_P_DF$theta_M, RE_theta_M_P_DF $mortality_P), 
                                   drop = TRUE),function(x)
                                           x[which.max(x$RE),]))

theta_M_contribution_to_RE_GG <- ggplot(RE_THETA_M_MAX_RE, 
       aes(x = theta_M/0.7, y= as.factor(mortality_P), fill = MtoH/RE)) + geom_raster() + 
        xlab(expression("Multiplier of transmission probability " * "(" *theta[M]* ")")) + 
        ylab("Disturbance intensity\nof primary vector") + 
        scale_fill_viridis(option = 'rocket', name = "Secondary \ncontribution to RE") + 
        scale_y_discrete(expand= c(0,0)) + 
        scale_x_continuous(expand = c(0,0)) + 
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 15))


theta_M_contribution_to_RE_GG / f_M_contribution_to_RE_GG + 
 plot_layout(guides = 'collect')

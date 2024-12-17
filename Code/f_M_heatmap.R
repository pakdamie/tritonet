eq_initial_list <- calculate_predisturb_initial_list()

# Running models and cleaning it up ---------------------------------------
modifier_values <- seq(0.05, 2.5, 0.10)
modified_f_M <- as.numeric(modifier * get_parameters("no_diff")[["f_M"]])
Mortality_P <- seq(0, 1, 0.05)
Mortality_P[1] <- 0.01
f_M_disturbance <- expand.grid(f_M = modified_f_M, 
                                   mortality_P = Mortality_P)


Mod_MortP_f_M <- Simulate_Model_Output(
  get_parameters("post_disturb"), 
  c("f_M", "mortality_P"), f_M_disturbance)


RE_f_M_P_List <- NULL


for (param in seq(1:nrow(f_M_disturbance))){
        
 param_changed <- get_parameters("post_disturb")
 param_changed["f_M"] <- f_M_disturbance[param,]$f_M
 param_changed["mortality_P"] <- f_M_disturbance[param,]$mortality_P   
        
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
        
  model_RE$f_M <- f_M_disturbance[param,]$f_M
  model_RE$mortality_P <- f_M_disturbance[param,]$mortality_P  
  
  RE_f_M_P_List[[param]] = model_RE
        
}

RE_f_M_P_DF <- do.call(rbind, RE_f_M_P_List)


RE_F_M_MAX_RE <- 
        do.call(rbind,
        lapply(split(RE_f_M_P_DF , 
        list(RE_f_M_P_DF$f_M,RE_f_M_P_DF $mortality_P), 
        drop = TRUE),function(x)
        x[which.max(x$RE),]))

f_M_contribution_to_RE_GG<- ggplot(RE_F_M_MAX_RE, 
       aes(x = f_M/0.02, y= as.factor(mortality_P), fill = MtoH/RE)) + geom_raster() + 
        xlab(expression("Multiplier of biting rate " * "(" *f[M]* ")")) + 
        ylab("Disturbance intensity\n of primary vector") + 
        scale_fill_viridis(option = 'rocket', name = "Secondary \ncontribution to RE") + 
        scale_y_discrete(expand= c(0,0)) + 
        scale_x_continuous(expand = c(0,0)) + 
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 15))

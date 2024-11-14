param_standard <-  
        c(b_H = 1/(1000), ##Human mortality rate
          b_P = 0.1, # P. Vector birth rate
          b_M = 0.1,  # S. Vector birth rate
          mu_H =1/(1000),  ##Human death rate
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
          disturbance_time = 8000,
          delta_T = 1,
          d = 0.05,
          prop = 1,
          mortality_P = 0.25,
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
Initial_List <- create_initial_states(param_standard, patch_num = 100)


adjacency_matrix <- as_adjacency_matrix(
        low_connectance_network,
        type = "both",
        attr = "weight",
        names = TRUE,
        sparse = FALSE
)
summed_prob <- rowSums(adjacency_matrix)
adjacency_matrix_adj <- sweep(adjacency_matrix, MARGIN = 1, summed_prob, `/`)


###Simulate the model with varying primary-removal 
model_output_mort_P <- NULL
for (p in seq(1,length(parameter_mortality_P_list))){
        
        model_output_mort_P[[p]] <- 
                discrete_trito_model_rcpp(
                        HS = Initial_List[[1]],
                        HI = Initial_List[[2]],
                        HR = Initial_List[[3]], 
                        PS = Initial_List[[4]],
                        PI = Initial_List[[5]],
                        MS = Initial_List[[6]],
                        MI = Initial_List[[7]],
                        adj = adjacency_matrix_adj,
                        param = parameter_mortality_P_list[[p]])    
        
}

## Calculate the RE

RE_multipatch_List <- NULL
RE_multipatch_patch_List <- NULL
for (p in seq(1,length(parameter_mortality_P_list))){
        
        RE_multipatch_List[[p]] <- cbind.data.frame(
          calculate_R_effective_discrete_net(
                (parameters = parameter_mortality_P_list[[p]]),
                 lists =model_output_mort_P[[p]],
                  8000),
                primary_removal = mortality_P[[p]])
        
        RE_multipatch_patch_List[[p]] <- cbind.data.frame(
                calculate_R_effective_discrete_patch(
                        (parameters = parameter_mortality_P_list[[p]]),
                        lists =model_output_mort_P[[p]],
                        8000),
                primary_removal = mortality_P[[p]])
}

RE_multipatch_DF <- do.call(rbind,RE_multipatch_List )
RE_multipatch_patch_DF <- do.call(rbind, RE_multipatch_patch_List)

RE_multipatch_patch_CV <- do.call(rbind.data.frame,lapply(split(RE_multipatch_patch_DF,list(RE_multipatch_patch_DF$time,
                                      RE_multipatch_patch_DF$primary_removal)),
        function(x) cbind(CV =sd(x$RE)/mean(x$RE),time = unique(x$time),
                         primary_removal= unique(x$primary_removal))))

ggplot(RE_multipatch_patch_CV ,aes( x= time - 8000, y = CV, 
                                    color = as.factor(primary_removal), 
                                    group = as.factor(primary_removal)))+ 
        geom_textline(aes(label =primary_removal )) +
        scale_color_viridis(discrete = TRUE, option = 'turbo') + 
        ylab("CV (Patch RE)") + xlab("Time") + 
        theme_classic() + 
        theme(legend.position = "none",
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 15))





ggplot(RE_multipatch_DF, 
       aes (x = time - 8000, y= RE,
            color = as.factor(primary_removal),
            group = as.factor(primary_removal))) + 
        geom_textline(aes(label =primary_removal )) +
        scale_color_viridis(discrete = TRUE, option = 'turbo') + 
        ylab("Network RE") + xlab("Time") + 
        theme_classic() + 
        theme(legend.position = "none",
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 15))

###The user inputs an FVec assuming that there is one patch (the function should
### change it )


adjacency_matrix <- matrix(c(0,1,1,
                             1,0,0,
                             1,0,0),nrow = 3, ncol = 3)

inf_states <- c("HI","PI","SI")
other_states <-  c("HS","HR", "PS", "SS","NH", "NP", "NS")

params <- data.frame(phi_P = 0.132, #transmission probability of p. vector
                     phi_S = 0.132 * 0.55, #transmission probability of s. vector
                     phi_H  = 0.116, 
                     a_P = 7, #daily biting rate of the p. vector
                     a_S = 7 * (0.55), 
                     gamma = 1/56,
                     mu_H = 1/27375, #Human death rate
                     mu_P = 0.05, #P. vector death rate
                     mu_S = 0.05, #S. vector death rate
                     amax = 2 * 10^-2 ,
                     k = 0.02,
                     am = 300,
                     cPS = 5e-4, #competitition effect of p.vector on s.vector
                     cSP = 1e-6,
                     a = 1)



tmp

F_matrix_test <- as.matrix(c("phi_P * a_P * HS1 * (PI1/NP1) + phi_S * a_S * HS1 *(SI1/NS1)",
                           "phi_H * a_P * PS1 * (HI1/NH1)",
                           "phi_H * a_S * SS1 * (HI1/NH1)"),
                           ncol = 1)




V_matrix_test <- as.matrix(c("(gamma + mu_H) * HI1",
                             "(mu_P + (amax/(a+exp(-k*(PI1 - am))))* PI1) + (cSP * SI1 * PI1) - MAT",
                             "(mu_S + (amax/(a+exp(-k*(PI1 - am))))* PI1)  + (cPS * PI1 * SI1)  - MAT"),
                           ncol = 1)



Calculate_R_effective <- function(adjacency_matrix, 
                           inf_states, other_states,
                           FVec, VVec, deSolve_Dataframe,
                           parameter_list){
        
        #Calculate the total number of patches
        num_patch <- nrow(adjacency_matrix)
        
        #How many infected classes that we are interested per patch?
        infected_full_states <- paste0(inf_states, rep(1:num_patch, 
                                          each = length(inf_states)))
        
        #F-MATRIX
        #Account for all other patches
        full_F_matrix <- NULL
        
        #Replace the number by the increment (HS1 to HS2)
        for (i in seq(2, num_patch + 1)){
                full_F_matrix[[i-1]]  <- unlist(lapply(FVec, function(x) 
                          gsub("[0-9]", i-1, x)))
        }
        
        #V-matrix
        # Row by row, figure out which column (other patches) are connected
        #to the patch of interest
        connected_patches = apply(adjacency_matrix, 1, function(x) which(x != 0))
         
        #This data.frame will allow us to figure out which equations to paste
        #into the model! Damie needs to generalize this.
        
        adjacency_df = NULL
        for (patch_interest in seq(1, num_patch)){
                adjacency_df [[patch_interest]] = cbind.data.frame(patch = patch_interest, 
                        connected_to_P = paste0("PI", connected_patches [[patch_interest]]),
                        connected_to_S = paste0("SI", connected_patches [[patch_interest]]))
        }  
        adjacency_df <- do.call(rbind, adjacency_df)   
        
         full_V_matrix <- NULL
         
         # Using information from the adjacency df, paste the necessary terms
         # for this patch: Oof fix this for generalizability
         for (patch_int in seq(1,num_patch)){
                 
                 patch_specific <- subset(adjacency_df, adjacency_df$patch == patch_int)
                 
                 added_PI_dispersal_term <- unlist(lapply(patch_specific$connected_to_P,
                                            function(x) paste0("(amax/(a+exp(-k*(PI1 - am))) * ", x,")")))
                 
                 added_SI_dispersal_term <- unlist(lapply(patch_specific$connected_to_S,
                                             function(x) paste0("(amax/(a+exp(-k*(PI1 - am))) * ",x,")")))
                 
                 full_PI_dispersal_term <- paste(added_PI_dispersal_term, collapse = " + ")
                 full_SI_dispersal_term <- paste(added_SI_dispersal_term, collapse = " + ")
                 
                 
                 V_mat_expression1 <- gsub("[0-9]", i, V_matrix_test[[1]])
                 V_mat_expression2 <- gsub("[0-9]", i, V_matrix_test[[2]])
                 V_mat_expression3 <- gsub("[0-9]", i, V_matrix_test[[3]])
                 
                 V_mat_expression2 <-sub("MAT", 
                                         full_PI_dispersal_term, V_mat_expression2)
                 
                 V_mat_expression3 <-sub("MAT", 
                                         full_SI_dispersal_term, 
                                         V_mat_expression3 )
                 
                 full_V_matrix[[patch_int]] <- list(V_mat_expression1,
                                            V_mat_expression2,
                                            V_mat_expression3)
                 
         }      
        
        #These are the full F and V-matrices that are now 
        #ready for analysis
        full_F_matrix <- unlist(full_F_matrix)
        full_V_matrix <- unlist(full_V_matrix)

   
        jacobian_F <- matrix(0, nrow = length(infected_full_states),
                                ncol = length(infected_full_states))
        
        jacobian_V <- matrix(0, nrow = length(infected_full_states),
                                ncol = length(infected_full_states))
        
        #Symbolically get the Jacobian matrix
        for(jac_ele in seq(1,length(infected_full_states))){
                
                state_interest <- infected_full_states[[jac_ele]]
                
                jacobian_F[,jac_ele] <- unlist(lapply(full_F_matrix, 
                                            function(x) 
                                            Deriv(x, state_interest)))
                
                jacobian_V[,jac_ele] <- unlist(lapply(full_V_matrix,
                                                      function(x)
                                                     Deriv(x, state_interest)))
                print(jacobian_V)
                
        }
 
       # for(time in seq(1, 1:nrow(deSolve_Dataframe))){
                hey=data.frame(params, tryit[2,])

                
                
#        }
        

        Evaluated_V <-apply(jacobian_V, 1:2,function(x) with(hey, eval(parse(text = x))))
        

        Evaluated_F <- apply(jacobian_F, 1:2, function(x) with(hey, eval(parse(text = x))))
                     
        #print(Evalulated_V)
        
        Evaluated_FV1 <- Evaluated_F  %*% solve(Evaluated_V)        
            
        

        R_Effective <- max(Re(eigen(Evaluated_FV1)$values))
           
        return(      Evaluated_V  )
        
        
        
}

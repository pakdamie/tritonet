###The user inputs an FVec assuming that there is one patch (the function should
### change it )


adjacency_matrix <- matrix(c(0,1,1,
                             1,0,0,
                             1,0,0),nrow = 3, ncol = 3)


R_effective_calculator <- (adjacency_matrix, I_states, Full_states,
                           FVec, VVec, deSolve_Dataframe,
                           parameter_list){
        
        
        
        #Calculate the total number of patches
        num_patch <- nrow(adjacency_matrix)
        
        states <- c("HI","PI","SI")
        non_infected_states <- c("HS","HR", "PS", "SS")
        full_states <- c("NH", "NP", "NS")
        
        
        total_state_variables <- c(states, non_infected_states,full_states)
                
        infected_full_states <- paste0(states, rep(1:num_patch, 
                                          each = length(states)))
        
        
        full_F_matrix <- NULL
        
        for (i in seq(2,num_patch + 1)){
        full_F_matrix[[i-1]]  <- unlist(lapply(F_matrix_test, function(x) 
                        gsub("[0-9]", i-1, x)))
        }
        
        ###
        
            
         connected_to = apply(adjacency_matrix, 1, 
                                  function(x) which(x!=0))
         adjacency_df = NULL
         for (i in seq(1,num_patch)){
                 adjacency_df [[i]] = cbind.data.frame(patch = i, 
                                                       connected_to_P = 
                                                               paste0("PI",
                                                               connected_to[[i]]),
                                                       connected_to_S = 
                                                               paste0("SI",
                                                                connected_to[[i]]))
         }  
         adjacency_df <- do.call(rbind, adjacency_df)   
         
         full_V_matrix <- NULL
         
         for (i in seq(1,num_patch)){
                 
                 patch_specific <- subset(adjacency_df, adjacency_df$patch == i)
                 
                 added_PI_dispersal_term <- unlist(lapply(patch_specific$connected_to_P,
                       function(x) paste0("(a * ",x,")")))
                 
                 added_SI_dispersal_term <- unlist(lapply(patch_specific$connected_to_S,
                                                          function(x) paste0("(a * ",x,")")))
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
                 
                 full_V_matrix[[i]] <- list(   V_mat_expression1,
                                               V_mat_expression2,
                                               V_mat_expression3)
                 
         }      
                
        full_F_matrix <- unlist(full_F_matrix)
        full_V_matrix <- unlist(full_V_matrix)

        jacobian_F <- matrix(0, nrow = length(infected_full_states),
                                ncol = length(infected_full_states))
        
        jacobian_V <- matrix(0, nrow = length(infected_full_states),
                       ncol = length(infected_full_states))
        
        for(jac_ele in seq(1,length( infected_full_states ))){
                
                state_interest <- infected_full_states[[jac_ele]]
                
                jacobian_F[,jac_ele] <- unlist(lapply(full_F_matrix, 
                                            function(x) 
                                                    Deriv(x, state_interest)))
                
                jacobian_V[,jac_ele] <- unlist(lapply(full_V_matrix,
                                                      function(x)
                                                     Deriv(x, state_interest)))
                
        }
 
        with(params, eval(parse(text =jacobian_V)))
        Evaluated_V <-apply(jacobian_V, 1:2,function(x) with(params, eval(parse(text =x))))
        Evaluated_F <- apply(jacobian_F, 1:2, function(x) with(params, eval(parse(text =x))))
                     
        Evaluated_FV1 <- Evaluated_F  %*% solve(Evaluated_V)        
            
        R_Effective <- max(Re(eigen(  Evaluated_FV1 )$values))
           
        rowSums(Evaluated_FV1)      
        
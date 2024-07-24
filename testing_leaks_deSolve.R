###Testing function to insure that there are no leaks

testing_leaks_deSolve<- function(test,num_patches){
      
        num_patches = 50
        new_network_null <-  erdos.renyi.game(
                num_patches,
               0.25,
                type = c("gnp"),
                directed = FALSE,
                loops = FALSE
        )
   
        adj_matrix_null <- as_adjacency_matrix(new_network_null, sparse = FALSE)
        
        ###The initial conditions
        initial_y_null <- c(HS = rep(1,num_patches),
                       HI = rep(1,num_patches),
                       HR = rep(1,num_patches),
                       PS = rep(5,num_patches),
                       PI = rep(1,num_patches),
                       SS = rep(5,num_patches),
                       SI = rep(1,num_patches))
        
        
        initial_y <- c(HS = rep(1,num_patches),
                       HI = rep(1,num_patches),
                       HR = rep(1,num_patches),
                       PS = sample.int(num_patches),
                       PI = sample.int(num_patches),
                       SS =  sample.int(num_patches),
                       SI = sample.int(num_patches))
        
        
        ###Make all of the parameters 0
        parameters_null <- c(
                b_H = 0, #Human birth rate
                b_P = 0, #P.vector birth rate
                b_S = 0, #S. vector birth rate
                mu_H = 0, #Human death rate
                mu_P = 0, #P. vector death rate
                mu_S = 0, #S. vector death rate
                
                a_P = 0.5, #biting rate of the p. vector
                a_S = 0.5, #biting rate of the s.vector
                
                phi_P = 0.5, #transmission probability of p. vector
                phi_S = 0.5, #transmission probability of s. vector
                phi_H  = 0.5, #transmission probability of human
                
                # Recovery rate
                gamma = 0,  #recovery rate of infected human
                
                #competition coefficient
                c_PS = 0, #competitition effect of p.vector on s.vector
                c_SP = 0,  #competitition effect of s.vector on p.vector
        
                disp_max = 1
                )
        
        

        results_null <- data.frame(deSolve::lsoda(
                times = 1:100,
                y =  initial_y ,
                func = trito_metapop,
                parms = parameters_null,
                patch_num =50,
                disp_mat = adj_matrix_null 
                
        ))
        
        
        
     Number_of_Individuals = sum(initial_y)
     
     All_individuals = rowSums(results_null[,(2:ncol(results_null))])
      
     if(any(round(Number_of_Individuals,3) != round(All_individuals,3))== TRUE){
             print("LEAKING! WARNING")
             }else{
                             print("No Leaks! Good!")
                     }
             
     

     
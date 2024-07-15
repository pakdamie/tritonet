Simulator_function <- function(num_patches,
                               initial_values = "default",
                               parameter_values = "default",
                               times){
       
        new_network <-  erdos.renyi.game(
                num_patches,
                0.35,
                type = c("gnp"),
                directed = FALSE,
                loops = FALSE
        )
                new_network <- set_vertex_attr(new_network, "name", 
                                            value=seq(1,num_patches))
        
        
        
        
        adj_matrix <- as_adjacency_matrix(new_network, sparse = FALSE)

        initial_y <- c(HS = rep(10000,num_patches),
                       HI = rep(10,num_patches),
                       HR = rep(10000,num_patches),
                       PS = rep(100,num_patches),
                       PI = rep(100,num_patches),
                       SS = rep(100,num_patches),
                       SI = rep(100,num_patches))
        
        
        
        parameters_n <- c(
                b_H = 0.05, #Human birth rate
                b_P = 0.02, #P.vector birth rate
                b_S = 0.02, #S. vector birth rate
                mu_H = 0.02, #Human death rate
                mu_P = 0.4, #P. vector death rate
                mu_S = 0.4, #S. vector death rate
                
                a_P = 0.90, #biting rate of the p. vector
                a_S = 0.90, #biting rate of the s.vector
                
                phi_P = 0.90, #transmission probability of p. vector
                phi_S = 0.10, #transmission probability of s. vector
                phi_H  = 0.40, #transmission probability of human
                
                # Recovery rate
                gamma = 1/7,  #recovery rate of infected human
                
                #competition coefficient
                c_PS = 0.004, #competitition effect of p.vector on s.vector
                c_SP = 0.00001  #competitition effect of s.vector on p.vector
        )
        
        df.contact <- adj_matrix * 0.02
        
        results <- data.frame(deSolve::lsoda(
                times = 1:5,
                y = initial_y ,
                func = trito_metapop,
                parms = parameters_n,
                patch_num =  num_patches ,
                disp_mat =df.contact 
                
        ))
        
        primary_infected <- cbind(time=results[,'time'],
        results [ , grepl( "PI" , names( results  ) ) ])
        primary_infected_Melted <- melt(primary_infected, id.vars ='time')
        primary_infected_Melted$Patch_num <- parse_number(as.character(primary_infected_Melted$variable))
        
       
      df_network <- fortify(new_network)
      nodes_coord <-  distinct(df_network,x,y,name)
      
      full_PI <- left_join(nodes_coord, primary_infected_Melted, by=c("name" = "Patch_num"))
      

      a<- ggplot(full_PI, aes(x =x, y=y, color = value))+geom_point(size = 4)+transition_time(time)
}

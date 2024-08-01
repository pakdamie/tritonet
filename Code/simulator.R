Simulator_function <- function(num_patch,
                               connectance,
                               max_distance = 20,
                               coverage,
                               frequency,
                               initial_values = "default",
                               parameter_values = "default",
                               disturbance = 'default', 
                               times){
        
        adjacency_matrix <- simulate_final_adjacency_matrix(24601, num_patch,
                                        connectance,max_distance)
        
        g9 <- graph_from_adjacency_matrix(adjacency_matrix , weighted=TRUE,
                                          mode="plus", diag=FALSE)
        
        
        sames <- sample(seq(1,100),num_patch, replace = TRUE)
        
        
        ###The initial conditions
        initial_y <- c(HS = sample(seq(1,1000),num_patch, replace = TRUE),
                       HI = rep(5,num_patch),
                       HR = rep(0,num_patch),
                       PS =  sames ,
                       PI = rep(5,num_patch),
                       SS =  sames ,
                       SI = rep(5,num_patch))
        

        parameters_null <- c(
                b_H = 1/27375, #Human birth rate
                b_P = 0.05, #P.vector birth rate
                b_S = 0.05, #S. vector birth rate
                mu_H = 1/27375, #Human death rate
                mu_P = 0.01, #P. vector death rate
                mu_S = 0.01, #S. vector death rate
                
                a_P = 4, #biting rate of the p. vector
                a_S = 2, #biting rate of the s.vector
                
                phi_P = 0.00008, #transmission probability of p. vector
                phi_S = 0.00004, #transmission probability of s. vector
                phi_H  = 0.5, #transmission probability of human
                
                # Recovery rate
                gamma = 1/56,  #recovery rate of infected human
                
                #competition coefficient
                c_PS = 0.0001, #competitition effect of p.vector on s.vector
                c_SP = 0.0001,  #competitition effect of s.vector on p.vector
                
                a_max = 0.0001,
                k = 0.001,
                a_0 = 1000, 
                
                lambda = 3     )
        
        
        
       chosen_patch <-  matrix(sample(seq(1, num_patch),floor(num_patch * coverage),
               replace = FALSE))
        
        
        disturbance <- data.frame(var = c(chosen_compartments( chosen_patch )),
                       value = c(rep(0.55,2 * length(chosen_patch)),
                                 rep(0.80, 2* length(chosen_patch))),
                       method = c(rep( "mult", 4 * length(chosen_patch) )))

        
       a<-  disturbance[rep(seq_len(nrow(disturbance)), length(seq(1,100,frequency))), ]
        
        
        a$time= rep(seq(1,100,frequency) ,each = nrow(disturbance))

        results <- data.frame(ode(
                times = seq(1,100,1),
                y =   initial_y  ,
                func = model_ross_trito_metapopulation,
                parms = parameters_null,
                patch_num =  num_patch ,
                adj_matrix = adjacency_matrix,
                events = list(data = a),
                method = 'lsodar'
                
        ))
        
        
        
        
        
        ### ALL SUBSEQUENT CODE AFTER IS PLACED SOMEWHERE ELSE.
        All_individuals = rowSums(results[,(2:ncol(results))])
        
        
        PSV_df <- calculate_PSV_ratio(results)
        
        df_network <- ggnetwork::fortify(   g9 )
        nodes_coord <-  distinct(df_network,x,y,name)
        
        nodes_coord$name <- as.numeric(nodes_coord$name) 
        
        full_PI <- left_join(nodes_coord, PSV_df, by=c("name" = "patch_num"))
        
        a_gif <- ggplot(data = df_network)+
                geom_segment(aes( x=x, xend = xend, y = y, yend = yend))+
                geom_point(data=full_PI,aes(x=x,y=y,color =value),size = 10)+
                scale_color_viridis(option = 'turbo', name = "Dominance")+                theme_void()+
                labs(title = 'Values at {(as.integer(frame_time))}')+
                transition_time(time)
        
        animate(a_gif, height = 800, width =800)
        anim_save("Gapminder_example.gif")
}
        

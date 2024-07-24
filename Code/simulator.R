Simulator_function <- function(num_patches,
                               initial_values = "default",
                               parameter_values = "default",
                               disturbance = 'default', 
                               times){
        
        ###We create a scale-free network
        new_network <-  erdos.renyi.game(
                num_patches,
                0.1,
                type = c("gnp"),
                directed = FALSE,
                loops = FALSE
        )
        
        ###Giving a specific name to the patches name
        new_network <- set_vertex_attr(new_network, "name", 
                                       value=seq(1,num_patches))
        
        ###This gives me the adjacency matrix
        adj_matrix <- as_adjacency_matrix(new_network, sparse = FALSE)
        
        ###The initial conditions
        initial_y <- c(HS = rep(1,num_patches),
                       HI = rep(5,num_patches),
                       HR = rep(0,num_patches),
                       PS = sample.int(num_patches),
                       PI = sample.int(num_patches),
                       SS =  sample.int(num_patches),
                       SI = sample.int(num_patches))
        
       sum(    initial_y )
        
        parameters_null <- c(
                b_H = 0, #Human birth rate
                b_P = 0, #P.vector birth rate
                b_S = 0, #S. vector birth rate
                mu_H = 0, #Human death rate
                mu_P = 0, #P. vector death rate
                mu_S = 0, #S. vector death rate
                
                a_P = 0.0, #biting rate of the p. vector
                a_S = 0.0, #biting rate of the s.vector
                
                phi_P = 0.0, #transmission probability of p. vector
                phi_S = 0.0, #transmission probability of s. vector
                phi_H  = 0.0, #transmission probability of human
                
                # Recovery rate
                gamma = 0.0,  #recovery rate of infected human
                
                #competition coefficient
                c_PS = 0.0, #competitition effect of p.vector on s.vector
                c_SP = 0.0,  #competitition effect of s.vector on p.vector
                
                disp_max = 1
        )
        
        df.contact <- adj_matrix 
        
        results <- data.frame(ode(
                times = seq(1,200,1),
                y =   initial_y  ,
                func = trito_metapop,
                parms = parameters_null,
                patch_num =  num_patches ,
                disp_mat = df.contact 
                
        ))
        All_individuals = rowSums(results[,(2:ncol(results))])
        
        
        PSV_df <- calculate_PSV_ratio(results)
        df_network <- fortify(new_network)
        nodes_coord <-  distinct(df_network,x,y,name)
        
        full_PI <- left_join(nodes_coord, PSV_df, by=c("name" = "Patch_num"))
        
        a_gif <- ggplot(data = df_network)+
                geom_segment(aes( x=x, xend = xend, y = y, yend = yend))+
                geom_point(data=full_PI,aes(x=x,y=y,color = log(value),size = log(value)))+
                scale_color_viridis(option = 'turbo')+
                scale_size_continuous(range=c(0.05,15))+
                theme_void()+
                labs(title = 'Values at {(as.integer(frame_time))}')+
                theme(legend.position = 'none')+
                transition_time(time)
        
        animate(a_gif, height = 800, width =800,rewind = TRUE)
        anim_save("Gapminder_example.gif")
}
        
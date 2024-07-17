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
        
        initial_y <- c(HS = rep(100,num_patches),
                       HI = rep(0,num_patches),
                       HR = rep(0,num_patches),
                       PS = sample.int(num_patches , replace = TRUE) * 100,
                       PI = rep(0,num_patches),
                       SS = sample.int(num_patches , replace = TRUE) * 100,
                       SI = rep(0,num_patches))
        
        
        
        parameters_n <- c(
                b_H = 0.05, #Human birth rate
                b_P = 11, #P.vector birth rate
                b_S = 11, #S. vector birth rate
                mu_H = 1/80, #Human death rate
                mu_P = 0.02, #P. vector death rate
                mu_S = 0.04, #S. vector death rate
                
                a_P = 0, #biting rate of the p. vector
                a_S = 0, #biting rate of the s.vector
                
                phi_P = 0.40, #transmission probability of p. vector
                phi_S = 0.10, #transmission probability of s. vector
                phi_H  = 0.90, #transmission probability of human
                
                # Recovery rate
                gamma = 1/7,  #recovery rate of infected human
                
                #competition coefficient
                c_PS = 0, #competitition effect of p.vector on s.vector
                c_SP = 0  #competitition effect of s.vector on p.vector
        )
        
        df.contact <- adj_matrix * 0 
        
        results <- data.frame(deSolve::lsoda(
                times = 1:120,
                y = initial_y ,
                func = trito_metapop,
                parms = parameters_n,
                patch_num =  num_patches ,
                disp_mat =df.contact 
                
        ))
        
        primary_infected <- cbind(time=results[,'time'],
                                  results [ , grepl( "PS" , names( results  ) ) ])
        
        secondary_infected <- cbind(time= results[,'time'],
                                    results[,grepl("SS", names(results))])
        
        
        primary_mat <- as.matrix(primary_infected[,2:(num_patches+ 1)])
        
        
        secondary_mat <- as.matrix(secondary_infected[,2:(num_patches+ 1)])
        
        
        primary_secondary_ratio <- cbind.data.frame(time = results[,'time'],
                                                    primary_mat/secondary_mat)
        
        
        primary_secondary_ratio_Melted <- melt(primary_secondary_ratio, id.vars ='time')
        primary_secondary_ratio_Melted$Patch_num <- parse_number(as.character(primary_secondary_ratio_Melted$variable))
        
        
        
        
        df_network <- fortify(new_network)
        nodes_coord <-  distinct(df_network,x,y,name)
        
        full_PI <- left_join(nodes_coord, primary_secondary_ratio_Melted, by=c("name" = "Patch_num"))
        
        a_gif <- ggplot(data = df_network)+
                geom_segment(aes( x=x, xend = xend, y = y, yend = yend))+
                geom_point(data=full_PI,aes(x=x,y=y,color = value),size = 8)+
                scale_color_viridis(option = 'turbo')+
                
                theme_classic()+
                labs(title = 'Values at {(as.integer(frame_time))}')+
                transition_time(time)
        
        animate(a_gif, height = 800, width =800)
        anim_save("Gapminder_example.gif")
        
        
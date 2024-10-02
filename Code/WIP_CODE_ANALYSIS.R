### A WIP ANALYSIS SCRIPT THAT IS A PLACEHOLDER.
coverage = c(0.25,.50, 1)
frequency = c(5,25,50)



expandDF <- expand.grid(coverage = coverage,
                        frequency = frequency)


LOW_CONNECTANCE_LIST <-NULL
MED_CONNECTANCE_LIST <- NULL
HIGH_CONNECTANCE_LIST <- NULL

networks_interest_low <- replicate(1,simulate_spatial_network(sample(1:24601,1),20,0.05))
networks_interest_med<- replicate(100,simulate_spatial_network(sample(1:24601,1),20,0.10))
networks_interest_high<- replicate(100,simulate_spatial_network(sample(1:24601,1),20,0.20))

LOW_NETWORK_LIST <- NULL
MED_NETWORK_LIST <- NULL
HIGH_NETWORK_LIST <- NULL

for (network in seq(1,100)){
        
        for (i in seq(1,nrow(expandDF))){
                tmp <- 
                Simulate_model(
                user_network = networks_interest_low[[network]],
                max_distance = 20,
                coverage = expandDF$coverage[1],
                frequency = expandDF$frequency[1],
                species1 = 0.5, 
                species2 = 0.8,
                initial_values = "default",
                parameter_values = "default",
                disturbance = "yes",         
                disease_on = 'no',
                simulation_time = 10,
                intervention_time = 10,
                randomize = TRUE)
                
               MED_CONNECTANCE_LIST[[i]] <- 
                       Simulate_model(
                               user_network = networks_interest_med[[network]],
                               max_distance = 20,
                               coverage = expandDF$coverage[i],
                               frequency = expandDF$frequency[i],
                               species1 = 0.5, 
                               species2 = 0.8,
                               initial_values = "default",
                               parameter_values = "default",
                               disturbance = "yes",         
                               disease_on = 'no',
                               simulation_time = 365,
                               intervention_time = 365,
                               randomize = TRUE)
                
                HIGH_CONNECTANCE_LIST[[i]] <- 
                        Simulate_model(
                                user_network = networks_interest_high[[network]],
                                max_distance = 20,
                                coverage = expandDF$coverage[i],
                                frequency = expandDF$frequency[i],
                                species1 = 0.5, 
                                species2 = 0.8,
                                initial_values = "default",
                                parameter_values = "default",
                                disturbance = "yes",         
                                disease_on = 'no',
                                simulation_time = 365,
                                intervention_time = 365,
                                randomize = TRUE)
                
}
               
        
       LOW_NETWORK_LIST[[network]] <- LOW_CONNECTANCE_LIST 
       MED_NETWORK_LIST[[network]] <- MED_CONNECTANCE_LIST 
       HIGH_NETWORK_LIST[[network]] <- HIGH_CONNECTANCE_LIST 
}

plot_patchwork_abundance(LOW_CONNECTANCE_LIST)
plot_patchwork_abundance(result_list1=LOW_CONNECTANCE_LIST, result_list2 = HIGH_CONNECTANCE_LIST)


Simulated_dominance_region_L <- lapply(LOW_CONNECTANCE_LIST,function(x) 
        calculate_dominance_region(x[[1]]))
Simulated_dominance_region_M <- lapply(MED_CONNECTANCE_LIST,function(x) 
        calculate_dominance_region(x[[1]]))
Simulated_dominance_region_H <- lapply(HIGH_CONNECTANCE_LIST,function(x) 
        calculate_dominance_region(x[[1]]))


for(i in seq(1,length(Simulated_dominance_region_L ))){
        Simulated_dominance_region_H[[i]]$coverage = expandDF$coverage[i]
        Simulated_dominance_region_H[[i]]$frequency= expandDF$frequency[i]
        Simulated_dominance_region_H[[i]]$connectance= c("H")
        
        Simulated_dominance_region_M[[i]]$coverage = expandDF$coverage[i]
        Simulated_dominance_region_M[[i]]$frequency= expandDF$frequency[i]
        Simulated_dominance_region_M[[i]]$connectance= c("M")
        
        Simulated_dominance_region_L[[i]]$coverage = expandDF$coverage[i]
        Simulated_dominance_region_L[[i]]$frequency= expandDF$frequency[i]
        Simulated_dominance_region_L[[i]]$connectance= c("L")
}

Simulated_dominance_region_DF_H <- do.call(rbind, Simulated_dominance_region_H)
Simulated_dominance_region_DF_M <- do.call(rbind, Simulated_dominance_region_M)
Simulated_dominance_region_DF_L <- do.call(rbind, Simulated_dominance_region_L)

Simulated_dominance_region_DF_L$dominant <- ifelse(Simulated_dominance_region_DF_L$prop >0.5,
                                                   "P",
                                                   "S")
Simulated_dominance_region_DF_H$dominant <- ifelse(Simulated_dominance_region_DF_H$prop >0.5,
                                                   "P",
                                                   "S")
Simulated_dominance_region_DF_M$dominant <- ifelse(Simulated_dominance_region_DF_M$prop >0.5,
                                                   "P",
                                                   "S")


Simulated_dominance_region_ALL <- rbind(Simulated_dominance_region_DF_L,
                                        Simulated_dominance_region_DF_M,
                                        Simulated_dominance_region_DF_H)


ggplot(Simulated_dominance_region_DF_L ,
       aes(x = time, y = prop,color = dominant, group = 1)) + geom_line(size = 0.8)+
        facet_grid(coverage~frequency) + theme_classic()+
        geom_hline(yintercept = 0.5)+
        xlab("Time")+
        ylab("Primary species proportion")+
        theme(strip.background = element_blank())+
        ggtitle("LOW")
        
ggplot(Simulated_dominance_region_DF_M ,
       aes(x = time, y = prop,color = dominant, group = 1)) + geom_line(size = 0.8)+
        facet_grid(coverage~frequency) + theme_classic()+
        geom_hline(yintercept = 0.5)+
        xlab("Time")+
        ylab("Primary species proportion")+
        theme(strip.background = element_blank())+
        ggtitle("MED")

ggplot(Simulated_dominance_region_DF_H ,
       aes(x = time, y = prop,color = dominant, group = 1)) + geom_line(size = 1)+
        facet_grid(coverage~frequency) + theme_classic()+
        geom_hline(yintercept = 0.5)+
        xlab("Time")+
        ylab("Primary species proportion")+
        theme(strip.background = element_blank())+ggtitle("HIGH")


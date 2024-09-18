### A WIP ANALYSIS SCRIPT THAT IS A PLACEHOLDER.
coverage = c(0.10,0.25,.50,1)
frequency = c(5,25,50,100,200)

###
list_networks <- simulate_spatial_network(24601,30)


expandDF <- expand.grid(coverage = coverage,
                        frequency = frequency)


LOW_CONNECTANCE_LIST <-NULL
MED_CONNECTANCE_LIST <- NULL
HIGH_CONNECTANCE_LIST <- NULL
for (i in seq(1,nrow(expandDF))){
        
LOW_CONNECTANCE_LIST[[i]] <- 
        Simulate_model(
                user_network = list_networks[[1]],
                max_distance = 30,
                coverage = expandDF$coverage[i],
                frequency = expandDF$frequency[i],
                species1 = 0.6, 
                species2 = 0.9,
                initial_values = "default",
                parameter_values = "default",
                disturbance = "yes",         
                disease_on = 'no',
                simulation_time = 300,
                intervention_time = 300)
        
        

MED_CONNECTANCE_LIST[[i]] <- 
        Simulate_model(
                user_network = list_networks[[2]],
                max_distance = 30,
                coverage = expandDF$coverage[i],
                frequency = expandDF$frequency[i],
                species1 = 0.6, 
                species2 = 0.9,
                initial_values = "default",
                parameter_values = "default",
                disturbance = "yes",         
                disease_on = 'no',
                simulation_time = 300,
                intervention_time = 300)

HIGH_CONNECTANCE_LIST[[i]] <- 
        Simulate_model(
                user_network = list_networks[[11]],
                max_distance = 30,
                coverage = expandDF$coverage[i],
                frequency = expandDF$frequency[i],
                species1 = 0.6, 
                species2 = 0.9,
                initial_values = "default",
                parameter_values = "default",
                disturbance = "yes",         
                disease_on = 'no',
                simulation_time = 300,
                intervention_time = 300)

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


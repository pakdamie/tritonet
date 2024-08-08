###coverage versus frequency

coverage = seq(0.1,1,0.1)
frequency = seq(5,50,5)
connectance = c(0.02)
expandDF <- expand.grid(coverage = coverage,frequency =frequency,
                        connectance = connectance)


Simulated <- mcmapply(Simulator_function,
                              num_patch =100,
                      connectance = expandDF$connectance,
                      max_distance = 20,
                      coverage = expandDF$coverage,
                      frequency = expandDF$frequency,
                      species1 = 0.5, species2 = 0.75,
                      initial_values = "default",
                      parameter_values = "default",
                      disturbance = "yes",         
                      disease_on = 'no',
                      end_length = 100,
                      mc.cores = 3,
                      SIMPLIFY = FALSE
)




Simulated_dominance_region <- lapply(Simulated , function(x) calculate_dominance_region (x)[100,])

for(i in seq(1,length(Simulated_dominance_region ))){
        Simulated_dominance_region[[i]]$coverage = expandDF$coverage[i]
        Simulated_dominance_region[[i]]$frequency= expandDF$frequency[i]
        Simulated_dominance_region[[i]]$connectance= expandDF$connectance[i]
        
}

Simulated_dominance_region_DF <- do.call(rbind, Simulated_dominance_region )

ggplot(
        Simulated_dominance_region_DF , aes(x = coverage, y= frequency, fill = prop ))+
        geom_tile() + scale_fill_viridis(name = "Prop of patches\ndominated by p.vector", option= "viridis")+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0))+
        xlab("Coverage")+
        ylab("Frequency")+
        facet_wrap(~connectance)


lapply()

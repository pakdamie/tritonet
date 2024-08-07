plotter_dynamics <- function(results, sprayed_patches){

        human_all<- cbind(time=results[,'time'],
                                   results [ , grepl("HS", names(results)) |
                                                     grepl( "HI", names(results)) |
                                                     grepl( "HR", names(results))])
           
        
        human_infected<- cbind(time=results[,'time'],
                               results [ , grepl("HI", names(results))])
        
        vector_sus<- cbind(time=results[,'time'],
                          results [  grepl( "PS", names(results)) |
                                            grepl( "SS", names(results))])
        
        vector_inf <- cbind(time=results[,'time'],
                            results [  grepl( "PI", names(results)) |
                                               grepl( "SI", names(results))])
        
        vector_P <- cbind(time=results[,'time'],
                          results [  grepl( "PS", names(results)) |
                                             grepl( "PI", names(results))])
        
        vector_S <- cbind(time=results[,'time'],
                          results [  grepl( "SS", names(results)) |
                                             grepl( "SI", names(results))])
        

        melted_H_infected <- melt(human_infected,id.vars = 'time')
        melted_H_infected$sprayed <- "NO"
        melted_H_infected$sprayed[ melted_H_infected$variable %in% 
                                   paste0("HI",chosen_patch)] <- "YES"
        
        gg_hi<- ggplot(melted_H_infected, aes(x=time, y = value, color = variable,linetype = sprayed))+
                geom_line(size =1) + 
                scale_color_viridis(option = 'viridis',discrete = TRUE) +
                xlab("Time") + 
                ylab("Abundance")+
                ggtitle("Infected humans")+
                 theme_classic()+
                 theme(legend.position = 'none',
                      axis.text = element_text(size =14),
                      axis.title = element_text(size =15))
                
        
        melted_P_vec <- melt(    vector_P,id.vars = 'time')
        melted_P_vec$sprayed <- "NO"
        melted_P_vec$sprayed[ melted_P_vec$variable %in% 
                                      c(chosen_compartments( chosen_patch ))] <- "YES"
        
        gg_pv <- ggplot( subset(melted_P_vec,
                                melted_P_vec$sprayed == "YES"), aes(x=time, y =log10(value+1),color = variable))+
                geom_line() + 
                scale_color_viridis(option = 'viridis',discrete = TRUE) +
                xlab("Time") + 
                ylab("Abundance")+
                ggtitle("Primary vectors")+
                theme_classic()+
                ylim(0,3)+
                theme(legend.position = 'none',
                      axis.text = element_text(size =14),
                      axis.title = element_text(size =15))     
        
        

melted_S_vec <- melt(    vector_S,id.vars = 'time')
melted_S_vec$sprayed <- "NO"
melted_S_vec$sprayed[ melted_S_vec$variable %in% 
                              c(chosen_compartments( chosen_patch ))] <- "YES"

gg_sv <- ggplot( melted_S_vec, aes(x=time, y = log10(value+1), color = variable,linetype = sprayed))+
        geom_line() + 
        scale_color_viridis(option = 'viridis',discrete = TRUE) +
        xlab("Time") + 
        ylab("Abundance")+
        ggtitle("Secondary vectors")+
        theme_classic()+
        ylim(0,3)+
        theme(legend.position = 'none',
              axis.text = element_text(size =14),
              axis.title = element_text(size =15)) ;gg_sv


gg_pv + gg_sv
}

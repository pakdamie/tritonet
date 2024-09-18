###Proportion of time dominated by one species or not 

Calculate_time_dominance_PS <- function(dominance_region){
       splitted_list <-  split(Simulated_dominance_region_ALL,
                              ~coverage+frequency+connectance)
        
        prop_p <- do.call(rbind,
                          lapply(splitted_list, function(x) 
                                  (table(x$dominant)[2]/sum(table(x)))))
        variables_interest <- unique(do.call(rbind,
                                      lapply(splitted_list, function(x) 
                           x[,c("coverage", "frequency","connectance")])))
 
        
        return(cbind(variables_interest,prop_p))
}

Time_dominance_PS<- Calculate_time_dominance_PS(Simulated_dominance_region_ALL)

Time_dominance_PS$connectance= factor(Time_dominance_PS$connectance, levels=c('L',"M","H"))


ggplot(Time_dominance_PS,
       aes(x = as.factor(coverage), y= as.factor(frequency), fill = S))+
        geom_tile(color = 'black')+facet_wrap(~connectance)+
        xlab("Coverage") + 
        ylab("Frequency")+
        scale_fill_viridis(option = 'mako')+
        scale_x_discrete(expand = c(0,0))+
        scale_y_discrete(expand = c(0,0))+
        coord_equal()+
        theme_classic()+ 
        theme(strip.background = element_blank())




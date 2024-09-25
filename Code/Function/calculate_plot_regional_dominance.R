calculate_dominance_region <- function(results_ODE){
    
        calculated_PSV <- calculate_PSV_ratio(results_ODE)
        
        splitted_PSV <- split(calculated_PSV,calculated_PSV$time)
        
        primary_dom <- do.call(rbind, lapply( splitted_PSV, 
                                              function(x) nrow(x[x$value >= 1,])))
        secondary_dom <- do.call(rbind, lapply( splitted_PSV, 
                                              function(x) nrow(x[x$value < 1,])))
        
       return(data.frame(time = seq(1:length(primary_dom)),
                          prop = primary_dom/(primary_dom + secondary_dom)))
       
}




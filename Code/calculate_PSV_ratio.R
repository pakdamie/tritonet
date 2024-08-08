calculate_PSV_ratio <- function(results){
        
        results = results[[1]]
        num_patch <- (ncol(results)-1)/7
        primary_susceptible <- cbind(time=results[,'time'],
                                  results [ , grepl( "PS", names(results))])
        
        primary_infected <- cbind(time=results[,'time'],
                                  results [ , grepl( "PI", names(results))])
        
        secondary_susceptible <- cbind(time= results[,'time'],
                                    results[,grepl("SS", names(results))])
        
        secondary_infected <- cbind(time= results[,'time'],
                                    results[,grepl("SI", names(results))])
        
        
        primarys_mat <- as.matrix(primary_susceptible[,2:(num_patch+ 1)])
        primaryi_mat <- as.matrix(primary_infected[,2:(num_patch+ 1)])
        secondarys_mat <- as.matrix(secondary_susceptible[,2:(num_patch+ 1)])
        secondaryi_mat <- as.matrix(secondary_infected[,2:(num_patch+ 1)])
        
       
        primary_secondary_ratio =cbind.data.frame(time = results[,'time'], 
                                (primarys_mat + primaryi_mat )/
                                (secondarys_mat + secondaryi_mat ))
        
        primary_secondary_ratio_Melted <- melt(primary_secondary_ratio, id.vars ='time')
        
        primary_secondary_ratio_Melted$patch_num <- parse_number(as.character(primary_secondary_ratio_Melted$variable))
        
        return(primary_secondary_ratio_Melted)
}

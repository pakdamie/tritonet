#' Calculates the ratio of primary to secondary vectors within a patch over time
#' This function first calculates the primary to secondary vector ratio for each 
#' patch across time. This is done by summing up the entire primary vector
#' (both susceptible and infected) and dividing it by the secondary vector
#' (both susceptible and infected)
#'
#' @param results 
#'
#' @return A melted data.frame that should show the ratio of primary to secondary
#' vectors for each patch across time
#' @export
#'
#' @examples
#' 
#' 
#' 
calculate_PSV_ratio <- function(results){
        
        results = results[[1]] ###Get the first element out
        num_patch <- (ncol(results)-1)/7
        
        ###Pull out the necessary columns 
        primary_susceptible <- cbind(time=results[,'time'],
                                  results [ , grepl( "PS", names(results))])
        
        primary_infected <- cbind(time=results[,'time'],
                                  results [ , grepl( "PI", names(results))])
        
        secondary_susceptible <- cbind(time= results[,'time'],
                                    results[,grepl("SS", names(results))])
        
        secondary_infected <- cbind(time= results[,'time'],
                                    results[,grepl("SI", names(results))])
        
        
        primarys_mat <- as.matrix(primary_susceptible[,2:(num_patch + 1)])
        primaryi_mat <- as.matrix(primary_infected[,2:(num_patch + 1)])
        secondarys_mat <- as.matrix(secondary_susceptible[,2:(num_patch + 1)])
        secondaryi_mat <- as.matrix(secondary_infected[,2:(num_patch + 1)])
        
       
        ###This should be how you calculate the primary to secondary.
        primary_secondary_ratio =cbind.data.frame(time = results[,'time'], 
                                (primarys_mat + primaryi_mat )/
                                (secondarys_mat + secondaryi_mat ))
        
        primary_secondary_ratio_Melted <- melt(primary_secondary_ratio, 
                                               id.vars ='time')
        
        primary_secondary_ratio_Melted$patch_num <- parse_number(as.character(
                primary_secondary_ratio_Melted$variable))
        
        return(primary_secondary_ratio_Melted)
}

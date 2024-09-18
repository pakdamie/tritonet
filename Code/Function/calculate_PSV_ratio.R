#' MAYBE NOT FINAL? Calculates the ratio of primary to secondary vectors within a patch over time
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
        
        list_abundance <- calculate_vector_abundance_patch(results)
       
        ###This should be how you calculate the primary to secondary.
        primary_secondary_ratio <- cbind.data.frame(time = results[,1], 
                                (list_abundance[[1]] +  list_abundance[[2]])/
                                (list_abundance[[3]] +  list_abundance[[4]]))
        
        primary_secondary_ratio <- primary_secondary_ratio[,-2]
        
        colnames(primary_secondary_ratio)[2:ncol(primary_secondary_ratio)] <- 
                seq(2:ncol(primary_secondary_ratio))
        
        
        
        primary_secondary_ratio_Melted <- melt(primary_secondary_ratio, 
                                               id.vars ='time')
        
        primary_secondary_ratio_Melted$patch_num <- parse_number(as.character(
                primary_secondary_ratio_Melted$variable))
        
        return(primary_secondary_ratio_Melted)
}

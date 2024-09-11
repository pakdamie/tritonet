

#' Find the closest values
#'
#' @param vec1_want 
#' @param vec2_have 
#'
#' @return
#' @export
#'
#' @examples
get_closest_values_vecs <- function(vec1_want, vec2_have){
        index_closest_value <- NULL
        for (value in 1:length(vec1_want)){
        index_closest_value [[value]] <- which.min(abs(vec1_want[[value]] - vec2_have))
        
        }
        return(do.call(rbind, index_closest_value ))
        

}
        

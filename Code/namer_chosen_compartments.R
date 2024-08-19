#' Labels the group of interest and the patch of interest 
#'
#' @param chosen_patch The patches that are chosen to be targeted
#'
#' @return A vector that includes "PS1", "PI1", etc
#' @export
#'
#' @examples
namer_chosen_compartments <- function(chosen_patch){
        
     ###Primary susceptible
     PS   <- apply(chosen_patch,2,
                 function(x) paste0("PS", x), simplify = TRUE)
        
     PI   <- apply(chosen_patch,2,
                   function(x) paste0("PI", x), simplify = TRUE)  
     
     SS   <- apply(chosen_patch,2,
                   function(x) paste0("SS", x), simplify = TRUE)  
     SI   <- apply(chosen_patch,2,
                   function(x) paste0("SI", x), simplify = TRUE)  
     
     return(c(PS, PI, SS, SI))
}




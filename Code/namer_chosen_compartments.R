
###This functions give you the names of the compartment and patch 
### to be used for the event function in deSolve

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




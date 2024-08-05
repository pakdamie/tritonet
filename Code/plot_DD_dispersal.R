###The function that we can use to figure out
### how to the density-dependent dispersal function would look like
### a0 is the maximum dispersal rate (daily),
### k is the steepness parameter, 
### N0 is the vector abundance that lead to the 

plotter_DD_dispersal <- function(a_max,k,N0){
        N_V = seq(1,1e4)
        dd_result <- a_max/(1.00 +  exp(-k *(N_V - N0)))
        
        plot(N_V, dd_result)
        
}

plotter_DD_dispersal(1,1e-2,2500)

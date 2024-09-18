#' Plot the density-dependence dispersal function (miscellaneous, just for 
#' making schematic figure or curiosity's sake)
#'
#' @param a_max The maximum rate of dispersal 
#' @param k The steepness slope
#' @param N0 The abundance at which there is the inflection point
#'
#' @return A plot of the DD dispersal over the total vector population
#' @export
#'
#' @examples
plotter_DD_dispersal <- function(a_max,k,N0){
        N_V = seq(1,1e4)
        dd_result <- a_max/(1.00 +  exp(-k *(N_V - N0)))
        plot(N_V, dd_result)
        
}

#plotter_DD_dispersal(1,1e-2,2500)

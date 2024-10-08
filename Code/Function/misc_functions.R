#All the packages that need to be loaded for the project
load_packages <- function(){
library(dplyr)
library(igraph)
library(ggplot2)
library(gganimate)
library(reshape2) 
library(readr)
library(Rcpp)
library(viridis)
library(deSolve)
library(forcats)
library(patchwork)
library(parallel)
library(Matrix)
}

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

#' Labels both the compartment and the patch of interest
#'
#' @param chosen_patch The patches that are chosen to be targeted
#'
#' @return A vector that includes "PS1", "PI1", etc
#' @examples 
namer_chosen_compartments <- function(chosen_patch) {
        
        compartments <- c("PS", "PI", #primary susceptible and infected
                          "SS", "SI") #secondary susceptible and infected

        # Generate the names for each compartment
        result <- lapply(compartments, function(prefix) {
                apply(chosen_patch, 2, function(x) paste0(prefix, x), 
                      simplify = TRUE)
        })
        
        
        return(do.call(rbind,result))
}




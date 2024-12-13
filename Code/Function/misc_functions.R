# All the packages that need to be loaded for the project
load_packages <- function() {
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
  library(geomtextpath)
}; load_packages()

#' Get different parameters
#'
#' A collection of different parameters I might have used, "standard" is the
#' default and is the one I would want
#'
#' @param type 
#'
#' @return A data.frame of parameter values.
#' @export
#'
#' @examples
get_parameters <- function(type = "standard"){
  param_standard <- c(
    b_H = 1 / (1000), ## Human mortality rate
    b_P = 0.01, # P. Vector birth rate
    b_M = 0.01, # S. Vector birth rate
    mu_H = 1 / (1000), ## Human death rate
    f_P = 0.02, # Biting rate of the p. vector
    f_M = 0.020 * 0.75, # Biting rate of the s.vector
    theta_P = 0.70, # Transmission probability of p. vector
    theta_M = 0.70 * 0.75, # Transmission probability of s. vector
    theta_H = 0.50, # Transmission probability of human
    gamma = 1 / 90, # Recovery rate of infected human
    c_PM = 4e-6, ## Competition effect of p.vector on s.vector
    c_MP = 2e-6, ## Competition effect of s.vector on p.vector
    c_PP = 4.5e-6, ## Competition effect of p.vector on s.vector
    c_MM = 3e-6, ## Competition effect of s.vector on s.vector
    ntime = 365 * 50,
    disturbance_time = 365 * 25,
    delta_T = 1,
    prop = 1,
    mortality_P = 0.25, # This will change
    mortality_M = 1)
  
  param_nodiff <- c(
    b_H = 1 / (1000), ## Human mortality rate
    b_P = 0.01, # P. Vector birth rate
    b_M = 0.01, # S. Vector birth rate
    mu_H = 1 / (1000), ## Human death rate
    f_P = 0.02, # Biting rate of the p. vector
    f_M = 0.020, # Biting rate of the s.vector
    theta_P = 0.70, # Transmission probability of p. vector
    theta_M = 0.70, # Transmission probability of s. vector
    theta_H = 0.50, # Transmission probability of human
    gamma = 1 / 90, # Recovery rate of infected human
    c_PM = 4e-6, ## Competition effect of p.vector on s.vector
    c_MP = 2e-6, ## Competition effect of s.vector on p.vector
    c_PP = 4.5e-6, ## Competition effect of p.vector on s.vector
    c_MM = 3e-6, ## Competition effect of s.vector on s.vector
    ntime = 365 * 50,
    disturbance_time = 365 * 25,
    delta_T = 1,
    prop = 1,
    mortality_P = 0.25, # This will change
    mortality_M = 1)
  
  param_no_disturb <- c(
    b_H = 1 / (1000), ## Human mortality rate
    b_P = 0.01, # P. Vector birth rate
    b_M = 0.01, # S. Vector birth rate
    mu_H = 1 / (1000), ## Human death rate
    f_P = 0.020, # Biting rate of the p. vector
    f_M = 0.020 * 0.75, # Biting rate of the s.vector
    theta_P = 0.70, # Transmission probability of p. vector
    theta_M = 0.70 * 0.75, # Transmission probability of s. vector
    theta_H = 0.50, # Transmission probability of human
    gamma = 1 / 90, # Recovery rate of infected human
    c_PM = 4e-6, ## Competition effect of p.vector on s.vector
    c_MP = 2e-6, ## Competition effect of s.vector on p.vector
    c_PP = 4.5e-6, ## Competition effect of p.vector on s.vector
    c_MM = 3e-6, ## Competition effect of s.vector on s.vector
    ntime = ((365 * 25)) ,
    disturbance_time = 1e30,
    delta_T = 1,
    prop = 1,
    mortality_P = 0.25, # This will change
    mortality_M = 1)
  
 param_post_distub <-  c(
     b_H = 1 / (1000), ## Human mortality rate
     b_P = 0.01, # P. Vector birth rate
     b_M = 0.01, # S. Vector birth rate
     mu_H = 1 / (1000), ## Human death rate
     f_P = 0.02, # Biting rate of the p. vector
     f_M = 0.02 , # Biting rate of the s.vector
     theta_P = 0.70, # Transmission probability of p. vector
     theta_M = 0.70, # Transmission probability of s. vector
     theta_H = 0.50, # Transmission probability of human
     gamma = 1 / 90, # Recovery rate of infected human
     c_PM = 4e-6, ## Competition effect of p.vector on s.vector
     c_MP = 2e-6, ## Competition effect of s.vector on p.vector
     c_PP = 4.5e-6, ## Competition effect of p.vector on s.vector
     c_MM = 3e-6, ## Competition effect of s.vector on s.vector
     ntime = (365 * 10) ,
     disturbance_time = 1,
     delta_T = 1,
     prop = 1,
     mortality_P = 0.25, # This will change
     mortality_M = 1
   )
 
  
  if(type == "standard"){
    return(param_standard )
  }
  else if(type == "no_disturb"){
    return(param_no_disturb)
  }
  else if(type == "post_disturb"){
   return(param_post_distub)
  }
 else if(type =="no_diff"){
   return(param_nodiff)
 }
}


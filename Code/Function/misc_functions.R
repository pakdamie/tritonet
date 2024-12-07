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
}
load_packages()

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
plotter_DD_dispersal <- function(a_max, k, N0) {
  N_V <- seq(1, 1e4)
  dd_result <- a_max / (1.00 + exp(-k * (N_V - N0)))
  plot(N_V, dd_result)
}


#' Labels both the compartment and the patch of interest
#'
#' @param chosen_patch The patches that are chosen to be targeted
#'
#' @return A vector that includes "PS1", "PI1", etc
#' @examples
namer_chosen_compartments <- function(chosen_patch) {
  compartments <- c(
    "PS", "PI", # primary susceptible and infected
    "SS", "SI"
  ) # secondary susceptible and infected

  # Generate the names for each compartment
  result <- lapply(compartments, function(prefix) {
    apply(chosen_patch, 2, function(x) paste0(prefix, x),
      simplify = TRUE
    )
  })


  return(do.call(rbind, result))
}

#' Find the closest values
#'
#' @param vec1_want
#' @param vec2_have
#'
#' @return
#' @examples
get_closest_values_vecs <- function(vec1_want, vec2_have) {
  index_closest_value <- NULL

  for (value in 1:length(vec1_want)) {
    index_closest_value[[value]] <- which.min(
      abs(vec1_want[[value]] - vec2_have)
    )
  }
  return(do.call(rbind, index_closest_value))
}

#' Plot the spatial netowrk graph
#'
#' @param seed The random seed generator
#' @param num_patch  Number of patches that you want to simulate
#' @param connectance The connectance of the network
#' @param max_distance The maximum distance
#'
#' @return
#' @export
#'
#' @examples
plot_networkgraph <- function(seed, num_patch, connectance, max_distance) {
  adjacency_matrix <- simulate_final_adjacency_matrix(
    seed, num_patch, connectance, max_distance
  )

  g9 <- graph_from_adjacency_matrix(adjacency_matrix,
    weighted = TRUE,
    mode = "plus", diag = FALSE
  )

  plot(g9)
}


#' Calculate the competition effects of P. and S. vectors
#'
#' @param list Model output
#' @param param parameter list
#'
#' @return A data.frame of the effects of the competition
#' @export
#'
#' @examples
calculate_competition_effects <- function(list, param) {
  c_MP <- param["c_MP"]
  c_PM <- param["c_PM"]
  c_MM <- param["c_MM"]
  c_PP <- param["c_PP"]

  PS <- list[[4]] # Susceptible Primary
  PI <- list[[5]] # Infectious Primary
  MS <- list[[6]] # Susceptible Secondary
  MI <- list[[7]] # Infectious Secondary

  NP <- PS + PI
  NM <- MS + MI

  comp_P_intra <- c_PP * NP
  comp_P_inter <- c_MP * NM
  comp_M_intra <- c_MM * NM
  comp_M_inter <- c_PM * NP

  return(cbind.data.frame(
    comp_P_intra, comp_P_inter,
    comp_M_intra, comp_M_inter,
    PS, PI, NP,
    MS, MI, NM
  ))
}

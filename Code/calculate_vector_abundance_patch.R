#' Extract the vector columns from a desolve output
#'
#' @param desolve_df 
#'
#' @return List with primary (susceptible and infected) as well as secondary 
#' (susceptible and infected)
#' @export
#'
#' @examples
calculate_vector_abundance_patch <- function(desolve_df) {
 
  num_patch <- (ncol(desolve_df) - 1) / 7 # the number of patches should simply be
  # the number of columns in desolve_df (-1, to get rid of the time column) divided
  ### by the number of different stages (H, P, and S)

  ### pick out the columns
  primary_susceptible <- cbind(
    time = desolve_df[, "time"],
    desolve_df[, grepl("PS", names(desolve_df))]
  )

  primary_infected <- cbind(
    time = desolve_df[, "time"],
    desolve_df[, grepl("PI", names(desolve_df))]
  )

  secondary_susceptible <- cbind(
    time = desolve_df[, "time"],
    desolve_df[, grepl("SS", names(desolve_df))]
  )

  secondary_infected <- cbind(
    time = desolve_df[, "time"],
    desolve_df[, grepl("SI", names(desolve_df))]
  )


  primarys_mat <- as.matrix(primary_susceptible[, 1:(num_patch + 1)])
  primaryi_mat <- as.matrix(primary_infected[, 1:(num_patch + 1)])
  secondarys_mat <- as.matrix(secondary_susceptible[, 1:(num_patch + 1)])
  secondaryi_mat <- as.matrix(secondary_infected[, 1:(num_patch + 1)])

  list_vectors <- list(
    primarys_mat,
    primaryi_mat,
    secondarys_mat,
    secondaryi_mat
  )
  names(list_vectors) <-c("PS","PI","SS","SI")
  
  return(list_vectors)
}

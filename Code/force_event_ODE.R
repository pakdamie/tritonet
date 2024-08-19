#' Force a disturbance (used in conjuction for ODE)
#'
#' @param num_patch Number of patches that you want to simulate
#' @param frequency How many times are the patches disturbed
#' @param coverage The proportion (0-1) of the patches that is targeted
#' @param species1  How much of the primary vector species survive
#' @param species2 How much of the secondary vector species survive
#' @param end_time total length of the simulation
#'
#' @return
#' @export
#'
#' @examples
force_event_ODE <- function(num_patch, frequency, coverage,
                            species1, species2, end_time) {
        
  # Sample from all the different patches. 
  # If patch is 
  chosen_patch <- matrix(sample(seq(1, num_patch), 
                                floor(num_patch * coverage), replace = FALSE
  ))

  disturbance_df <- data.frame(
    var = c(namer_chosen_compartments(chosen_patch)),
    value = c(
      rep(species1, 2 * length(chosen_patch)),
      rep(species2, 2 * length(chosen_patch))
    ),
    method = c(rep("mult", 4 * length(chosen_patch)))
  )


  event_df <- disturbance_df[rep(seq_len(nrow( disturbance_df )),
                              length(seq(1,end_time,frequency))), ]

  event_df$time <- rep(seq(1, end_time, frequency), each = nrow(disturbance_df))
  return(event_df)
}

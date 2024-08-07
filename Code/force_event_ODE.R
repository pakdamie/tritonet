### This is the disturbance function that we incorporate into the
### desolve equation
### supply with num_patch which is the number of patches,
### frequency - how many times to disturb (day)
### coverage- the proportion of patches
### species1 - survivorship (How many should survive?)
### species2 - survivorship (How many should survive?)
### end_time - how long should the disturbance run for?


force_event_ODE <- function(num_patch, frequency, coverage,
                            species1, species2, end_time) {
  # Sample from the different number of patches
  chosen_patch <- matrix(sample(seq(1, num_patch), floor(num_patch * coverage),
    replace = FALSE
  ))

  # Sample from the different number of patches

  disturbance_df <- data.frame(
    var = c(namer_chosen_compartments(chosen_patch)),
    value = c(
      rep(species1, 2 * length(chosen_patch)),
      rep(species2, 2 * length(chosen_patch))
    ),
    method = c(rep("mult", 4 * length(chosen_patch)))
  )


  event_df <- disturbance[rep(seq_len(nrow(disturbance)),
                              length(seq(1,end_time,frequency))), ]


  event_df$time <- rep(seq(1, end_time, frequency), each = nrow(disturbance))
  return(event_df)
}

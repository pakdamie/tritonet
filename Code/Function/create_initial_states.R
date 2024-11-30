#' Create initial states to simulate the model
#'
#' Create the seven matrices for each of the seven compartments. The seven
#' matrices will have the number of rows that depend on time-steps and
#' the number of columns is the number of patches.
#'
#' @param param The parameters that you're choosing
#' @param patch_num The number of patches
#'
#' @return A list of the necessary matrices
#' @export
#'
#' @examples create_initial_states(param_standard, 100)

create_initial_states <- function(
    param,
    patch_num,
    initial_human = 1000,
    initial_vector_s = 1000,
    initial_vector_i = 10) {
  
  compartment_labels <-
    c(
      "HS_mat", "HI_mat", "HR_mat", # Humans (sus/inf/rec)
      "PS_mat", "PI_mat", # Primary (sus/inf)
      "MS_mat", "MI_mat" # Secondary(sus/inf)
    )

  # Using the above compartment label, we then create an empty matrix where
  # the row are the time-steps and the columns are the individual patches

  for (i in 1:length(compartment_labels)) {
    assign(
      compartment_labels[i],
      matrix(0,
        nrow = param["ntime"],
        ncol = patch_num
      )
    )
  }

  # Initial conditions
  
  # Humans
  HS_mat[1, ] <- rep(initial_human, patch_num)
  HI_mat[1, ] <- rep(1, patch_num)
  HR_mat[1, ] <- rep(0, patch_num)

  # Primary vectors
  PS_mat[1, ] <- rep(initial_vector_s, patch_num)
  PI_mat[1, ] <- rep(initial_vector_i, patch_num)

  # Secondary vectors
  MS_mat[1, ] <- rep(initial_vector_s, patch_num)
  MI_mat[1, ] <- rep(initial_vector_i, patch_num)

  return(list(
    HS_mat, HI_mat, HR_mat,
    PS_mat, PI_mat,
    MS_mat, MI_mat
  ))
}

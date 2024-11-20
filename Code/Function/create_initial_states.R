#' Create initial states to simulate the model
#'
#' @param param The parameters that you're choosing
#' @param patch_num The number of patches
#'
#' @return A list of the necessary matrices
#' @export
#'
#' @examples
create_initial_states <- function(param, patch_num){

  compartment_label <- c(
    "HS_mat", "HI_mat", "HR_mat",  # Humans (sus/inf/rec)
    "PS_mat", "PI_mat",  # Primary (sus/inf)
    "MS_mat", "MI_mat"   # Secondary(sus/inf)
 )

  # Using the above compartment label, we then create an empty matrix where 
  # the row are the time-steps and the columns are the individual patches

  for (i in 1:length(compartment_label)) {
    assign(compartment_label[i],
           matrix(0, nrow = 
             param['ntime'], ncol = patch_num)
          )
  }

  #Initial conditions (everyone starts out with the same number)
  HS_mat[1, ] <- rep(1000, patch_num) 
  HI_mat[1, ] <- rep(1, patch_num) 
  HR_mat[1, ] <- rep(0, patch_num)
  
  PS_mat[1, ] <- rep(1000, patch_num)
  PI_mat[1, ] <- rep(1, patch_num)
  MS_mat[1, ] <- rep(1000, patch_num)
  MI_mat[1, ] <- rep(1, patch_num)

  
  return(list(HS_mat, HI_mat, HR_mat, 
              PS_mat, PI_mat, 
              MS_mat, MI_mat))
}

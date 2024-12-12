#' Vary a parameter at a time for the parameter list
#'
#' Vary the variable of interest to investigate how the model would change.
#' For example, the most important variables to be examining is the effect
#' of the disturbance.
#'
#' @param parameter List of parameters used for the model
#' @param variable_interest Character: What parameter do you want to vary?
#' @param vector_value Data_frame: What are the different values of that variable you want 
#' to vary?
#'
#' @return A list of parameters to simulate model with.
#' @export
#'
#' @examples vary_one_parameter_value(parameter_standard, "mortality_P", seq(0,1,0.01))

vary_parameter_value <- function(parameter, variable_interest, vector_value){
  
  if(class(variable_interest) != "character") 
    stop("The variable_interest should be a data.frame of characters!")
  
  if(class(vector_value) != "data.frame") 
    stop("The vector_value should be a data.frame of numeric values!")
  
  parameter_list <- NULL      
  copied_param <- parameter
  
  #Replacing the old value with vector_value[val]

    for (val in 1:nrow(vector_value)){
      
      for (var in 1:length(variable_interest)){
        
      var_interest =  variable_interest[[var]]
      copied_param[var_interest] <- vector_value[val,var]
      parameter_list[[val]] <-  copied_param
    }
  }
  
  return(parameter_list)
}


#' Create initial states to simulate the model
#'
#' Create the seven matrices for each of the seven compartments. The seven
#' matrices will have the number of rows that depend on time-steps and
#' the number of columns is the number of patches.
#'
#' @param param The parameters that you're choosing
#' @param patch_num The number of patches
#' @return A list (with length 7) of the necessary matrices to simulate the model
#' @export
#'
#' @examples create_initial_states(param_standard, 100)

create_initial_states <- function(param, patch_num = 1, initial_human = 1000,
  initial_vector_s = 1000, initial_vector_i = 10) {
        
  compartment_labels <- c(
    "HS_mat", "HI_mat", "HR_mat", # Humans (sus/inf/rec)
    "PS_mat", "PI_mat", # Primary (sus/inf)
    "MS_mat", "MI_mat") # Secondary(sus/inf))
        
  # Using the above compartment label, we then create an empty matrix where
  # the row are the time-steps and the columns are the individual patches
        
  for (i in 1:length(compartment_labels)) {
    assign(compartment_labels[i], matrix(0,
      nrow = param["ntime"], ncol = patch_num))
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
    MS_mat, MI_mat)
   )
}


#' Simulate the model 
#'
#' The main function for simulating the outputs of the model.
#' @param parameter 
#' @param variable_interest 
#' @param vector_value 
#'
#' @return
#' @export
#'
#' @examples
Simulate_Model_Output <- function(parameter, variable_interest, vector_value){

  parameter_list <- vary_parameter_value(
    parameter, variable_interest, vector_value)
   
  model_output_list <- NULL
  for (param in seq(1:nrow(vector_value))){
          
  model_output_list[[param]] <- 
    discrete_trito_model_rcpp_ONEPATCH(
    HS = Initial_List[[1]],
    HI = Initial_List[[2]],
    HR = Initial_List[[3]],
    PS = Initial_List[[4]],
    PI = Initial_List[[5]],
    MS = Initial_List[[6]],
    MI = Initial_List[[7]],
    param = parameter_list[[param]])
  }
 return(model_output_list) 
}


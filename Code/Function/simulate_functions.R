#' Retrieve different data.frame of parameters for the model simulation.
#'
#' A collection of different parameters I use for analysis: "standard" is the
#' default.
#'
#' @param type The type of parameters to get. The default is "standard", but other
#' types could include: no_disturb, post_disturb, post_disturb, and no_diff.
#' @return (data.frame) A data.frame of parameter values to use for the model
#' @export
#'
#' @examples get_parameters("standard")
get_parameters <- function(type = "standard") {
  
  if (!(type %in% c(
    "standard",
    "no_diff",
    "worse_m",
    "nonesec",
    "better_m"
  ))) 
    {
    stop("Ensure that you're using a valid type")
  }
  ### The standard parameter that is in the supplement
  param_standard <- c(
    b_H = 1 / (365 * 70), # Human mortality rate
    b_P = 1 / 10, # P. Vector birth rate
    b_M = 1 / 10, # S. Vector birth rate
    mu_H = 1 / (365 * 70), # Human death rate
    f_P = 0.25, # Biting rate of the p. vector
    f_M = 0.25, # Biting rate of the s.vector
    theta_P = 0.25, # Transmission probability of p. vector
    theta_M = 0.25 * 0.80, # Transmission probability of s. vector
    theta_H = 0.25, # Transmission probability of human
    gamma = 1 / 90, # Recovery rate of infected human
    c_PM = 3e-4, # Competition effect of p.vector on s.vector
    c_MP = 3e-6, # Competition effect of s.vector on p.vector
    c_PP = 4.5e-4, # Competition effect of p.vector on s.vector
    c_MM = 2.5e-4, # Competition effect of s.vector on s.vector
    mu_V = 1 / 20, # Vector mortality rate
    ntime = (365 * 100), # Number of times
    disturbance_time = (365 * 50),
    infection_time = 1e30,
    delta_T = 1,
    prop = 1,
    mortality_P = 0.50, # This will change
    mortality_M = 1
  )

  ### Secondary vectors exist (and compete), but does not contribute
  # to the RE in anyway
  param_nonesec <- param_standard
  param_nonesec[c("f_M", "theta_M")] <- 0
 
   ### Secondary vectors are even worse than the standard
  param_worse_m <- param_standard
  param_worse_m["f_M"] <- param_standard["f_M"] * 0.50
  param_worse_m["theta_M"] <- param_standard["theta_M"] * 0.50
  
  ### Secondary vectors have the same competence in transmission, but
  ## still worse competitor
  param_nodiff <- param_standard
  param_nodiff[c("f_M")] <- param_standard[c("f_P")]
  param_nodiff[c("theta_M")] <- param_standard[c("theta_P")]

  ## Secondary vectors are better transmitter but nothing is changed
  ## about their competition
  param_better_m <- param_standard
  param_better_m["f_M"] <- param_standard["f_M"] * 1.25
  param_better_m["theta_M"] <- param_standard["theta_M"] * 1.25

  # Depending on type
  param <- switch(type,
    "standard" = param_standard,
    "post_disturb" = param_post_disturb,
    "no_diff" = param_nodiff,
    "better_m" = param_better_m,
    "worse_m" = param_worse_m,
    "nonesec" = param_nonesec
  )

  return(param)
}


#' Vary a parameter at a time for the parameter list
#'
#' Vary the variable of interest to investigate how the model would change.
#' For example, the most important variables to be examining is the effect
#' of the disturbance. You assign it to the `variable_interest` and give
#' the different values that you are interested in. 
#'
#' @param parameter (list) of parameters used for the model
#' @param variable_interest (character): What parameter do you want to vary?
#' @param vector_value (data.frame): What are the different values of that variable you want to vary?

#' @return A list of parameters to simulate model with.
#' @export
#'
#' @examples vary_one_parameter_value(parameter_standard, "mortality_P", seq(0, 1, 0.01))
vary_parameter_value <- function(parameter,
                                 variable_interest,
                                 vector_value) {
  # Error messages if input is incorrect
  if (class(variable_interest) != "character") {
    stop("The variable_interest should be a data.frame of characters!")
  }
  if (class(vector_value) != "data.frame") {
    stop("The vector_value should be a data.frame of numeric values!") 
  }

  # Initialize an empty list to store modified parameter data.frames
  parameter_list <- NULL

  # Create a copy of the parameter object to modify later
  copied_param <- parameter

  # Iterate over each row in vector_value to update the parameter values
  for (val in 1:nrow(vector_value)) {
    # Iterate over each variable in variable_interest
    for (var in 1:length(variable_interest)) {
      var_interest <- variable_interest[[var]] # Extract variable name
      copied_param[var_interest] <- vector_value[val, var] # Assign new value
      parameter_list[[val]] <- copied_param # Store modified parameter set in the list
    }
  }

  # Return the list containing modified parameter sets
  return(parameter_list)
}


#' Create initial states to simulate the model
#'
#' Create the seven matrices for each of the seven compartments. The seven
#' matrices will have the number of rows that depend on time-steps and
#' the number of columns is the number of patches.
#' 
#' By default, we do not have infection started ("No"), but you can convert it
#' to ("Yes") to have infection started prior to the perturbation (make it
#' be endemic)
#'
#' @param param (data.frame) The parameters that you're choosing
#' @param infection_start (character) Do you initiate infection prior to disturbance? (endemic)
#' @param initial_human (numeric) How many susceptible humans do you start out with
#' @param initial_vector_s (numeric) How many susceptible vectors do you start out with
#' @param initial_vector_i (numeric) How many infected vectors do you start out with 
#' @return (list) With length of length 7 of the necessary matrices to simulate the model
#' @export
#'
#' @examples create_initial_states(param_standard, 100)
create_initial_states <- function(param,
                                  infection_start = "No",
                                  initial_human = 1000,
                                  initial_vector_s = 2490,
                                  initial_vector_i = 10) {
  compartment_labels <- c(
    "HS_mat", "HI_mat", "HR_mat", # Humans (Sus/Inf/Rec)
    "PS_mat", "PI_mat", # Primary (Sus/inf)
    "MS_mat", "MI_mat" # Secondary(Sus/inf))
  )

  # Using the above compartment label,
  # we then create an empty matrix where
  # the row are the time-steps and the columns are the individual patches

  for (i in 1:length(compartment_labels)) {
    assign(compartment_labels[i], matrix(0,
      nrow = param["ntime"], ncol = 1
    ))
  }

  # If we don't have infection ("No"), just start out with the first
  # values, if we have infection initiate infection in the
  # vectors- uncomment if you want to start just in human

  initial_vector_s <- ifelse(infection_start == "No", 2500, 2490)
  initial_vector_i <- ifelse(infection_start == "No", 0, 10)
  # initial_h_i <- ifelse(infection_start == "No", 0, 0)


  # Initial conditions
  # Humans
  HS_mat[1, ] <- 1000
  HI_mat[1, ] <- 0
  HR_mat[1, ] <- 0

  # Primary vectors
  PS_mat[1, ] <- initial_vector_s
  PI_mat[1, ] <- initial_vector_i

  # Secondary vectors
  MS_mat[1, ] <- initial_vector_s
  MI_mat[1, ] <- initial_vector_i

  # Return a list where there are 7 matrices

  return(list(
    HS_mat, HI_mat, HR_mat,
    PS_mat, PI_mat,
    MS_mat, MI_mat
  ))
}
 


#' Simulate the main model
#'
#' The main function for simulating the one-patch, two vector, one host model.
#' 
#' @param parameter (data.frame) Parameters of interest (retrieved from `get_parameters`)
#' @param variable_interest (character) Variable that you are interested in manipulating 
#' @param vector_value (data.frame) The values of the variable that you are interested in
#'
#' @return (list) A list with seven elements OR a list with `variable_interest` 
#' number of elements with each element having seven elements each. 
#' @export
#'
#' @examples Simulate_Model_Output(get_parameters("standard))
Simulate_Model_Output <- function(parameter, 
                                  infection_start = "No",
                                  variable_interest = NA, 
                                  vector_value = NA) {
  
  
  #Warn the user that there's infection or not
  print(paste0("You are assuming ", infection_start, " infection"), sep = " ")

  #Create the initial conditions
  Initial_List <- create_initial_states(
    parameter,
    infection_start = infection_start
  )

  #If there's going to be any varying of the parameter,
  #use `vary_parameter_value` to make a list of parameters
  #to run
  
  if (any(is.na(variable_interest)) == FALSE) {
    parameter_list <- vary_parameter_value(
      parameter, variable_interest, vector_value
    )

    model_output_list <- NULL

    for (param in seq(1:nrow(vector_value))) {
      model_output_list[[param]] <-
        model_vectors_host(
          HS = Initial_List[[1]],
          HI = Initial_List[[2]],
          HR = Initial_List[[3]],
          PS = Initial_List[[4]],
          PI = Initial_List[[5]],
          MS = Initial_List[[6]],
          MI = Initial_List[[7]],
          param = parameter_list[[param]]
        )
    }
  } else { #If the parameters are not being varied
    model_output_list <-
      model_vectors_host(
        HS = Initial_List[[1]],
        HI = Initial_List[[2]],
        HR = Initial_List[[3]],
        PS = Initial_List[[4]],
        PI = Initial_List[[5]],
        MS = Initial_List[[6]],
        MI = Initial_List[[7]],
        param = parameter
      )
  }

  return(model_output_list)
}

#' Retrieve different data.frame of parameters for the model simulation.
#'
#' A collection of different parameters I might have used, "standard" is the
#' default. Other analyses may require two set of parameters. 
#'
#' @param type The type of parameters to get. The default is "standard", but other
#' types could include: no_disturb, post_disturb, post_disturb, and no_diff.
#'
#' @return (data.frame) A data.frame of parameter values to.  
#' @export
#'
#' @examples
get_parameters <- function(type = "standard"){
  
  if(!(type %in% c("standard","no_disturb","post_disturb", "no_diff","worse_m","nonesec",
                   "better_m"))) 
    stop("The type should either be `standard`, `no_disturb`, `post_disturb`,
         `no_diff`")
  
  param_standard <- c(
    b_H = 1 /(365 * 70), ## Human mortality rate
    b_P = 0.01, # P. Vector birth rate
    b_M = 0.01, # S. Vector birth rate
    mu_H = 1 /(365 * 70), ## Human death rate
    f_P = 0.25, # Biting rate of the p. vector
    f_M = 0.25 * 0.75, # Biting rate of the s.vector
    theta_P = 0.20, # Transmission probability of p. vector
    theta_M = 0.20 * 0.75, # Transmission probability of s. vector
    theta_H = 0.50, # Transmission probability of human
    gamma = 1 / 100, # Recovery rate of infected human
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
  
  
  ###Secondary vectors exist (and compete), but does not contribute
  #to the RE in anyway 
  param_nonesec <- param_standard 
  param_nonesec[c("f_M", "theta_M")]<- 0
  
  
  ###Secondary vectors have the same competence in transmission, but
  ## still worse competitor
  param_nodiff <- param_standard 
  param_nodiff[c("f_M")]<- param_standard[c("f_P")]
  param_nodiff[c("theta_M")]<- param_standard[c("theta_P")]
  
  ###Secondary vectors are the worst 
  param_worse_m<- param_standard 
  param_worse_m["f_M"] <- param_standard["f_M"] * 0.10
  param_worse_m["theta_M"] <- param_standard["theta_M"] * 0.10
  
  ##Secondary vectors better but they are weaker competitors 
  param_better_m <-   param_standard 
  param_better_m["f_M"] <- param_standard["f_M"] * 1.10
  param_better_m["theta_M"] <- param_standard["theta_M"] * 1.10
  
  if(type == "standard"){
    return(param_standard)
  }
  else if(type == "no_disturb"){
    return(param_no_disturb)
  }
  else if(type == "no_diff"){
    return(param_nodiff)
  }
  else if (type == "better_m"){
    return( param_better_m)
  }
  else if (type == "worse_m"){
    return(param_worse_m)
  }
  else if(type == "nonesec"){
    return(param_nonesec)
  }
}

#' Vary a parameter at a time for the parameter list
#'
#' Vary the variable of interest to investigate how the model would change.
#' For example, the most important variables to be examining is the effect
#' of the disturbance.
#'
#' @param parameter List of parameters used for the model
#' @param variable_interest Character: What parameter do you want to vary?
#' @param vector_value Data_frame: What are the different values of that variable you want to vary?
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
create_initial_states <- function(param, 
                                  patch_num = 1, 
                                  initial_human = 1000, 
                                  initial_vector_s = 2490, 
                                  initial_vector_i = 10) {
        
  compartment_labels <- c(
    "HS_mat", "HI_mat", "HR_mat", # Humans (Sus/Inf/Rec)
    "PS_mat", "PI_mat", # Primary (Sus/inf)
    "MS_mat", "MI_mat") # Secondary(Sus/inf))
        
  # Using the above compartment label,
  # we then create an empty matrix where
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

#' Calculate the equilibrium states prior to the disturbance 
#'
#' @returns
#' @export
#'
#' @examples
simulate_predisturb_initial_list <- function(){
  
  param_no_disturb <- get_parameters("no_disturb")
  Initial_List<- create_initial_states(get_parameters("standard"))
  Mod_Predisturb <- discrete_trito_model_rcpp_ONEPATCH(
    HS = Initial_List[[1]],
    HI = Initial_List[[2]],
    HR = Initial_List[[3]],
    PS = Initial_List[[4]],
    PI = Initial_List[[5]],
    MS = Initial_List[[6]],
    MI = Initial_List[[7]],
    param = param_no_disturb)
  
  Predisturb_EqStates <- lapply(Mod_Predisturb,
                               function(x) x[((365 * 25)),])
  
  eq_initial_list<- create_initial_states(get_parameters("post_disturb"))
  eq_initial_list[[1]][1] <- Predisturb_EqStates[[1]]
  eq_initial_list[[2]][1] <- Predisturb_EqStates[[2]]
  eq_initial_list[[3]][1] <- Predisturb_EqStates[[3]]
  eq_initial_list[[4]][1] <- Predisturb_EqStates[[4]]
  eq_initial_list[[5]][1] <- Predisturb_EqStates[[5]]
  eq_initial_list[[6]][1] <- Predisturb_EqStates[[6]]
  eq_initial_list[[7]][1] <- Predisturb_EqStates[[7]]
  
  return(eq_initial_list)
}

#' Simulate the main model
#'
#' The main function for simulating the one-patch, two vector, one host model. 
#' @param parameter The data.frame of parameters (retrieved from `get_parameters`)
#' @param variable_interest The variables that you are interested in manipulating (character)
#' @param vector_value The data.frame of values that you are interested in 
#'
#' @return (list) A list with seven elements. Each element includes a matrix with the
#' rows being the time of simulation
#' @export
#'
#' @examples
Simulate_Model_Output <- function(parameter, variable_interest, vector_value){

  Initial_List <- create_initial_states(parameter)
  
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

#' Simulate with different parameters post-disturbance
#'
#' When we want to see how parameters change RE, it is important that
#' we are doing the correct comparison. For this, we initialize with the 
#' same equilibrium. 
#'
#' @param variable_interest 
#' @param vector_value 

Simulate_Model_Output_PostD <- function(variable_interest, vector_value){
  
  param_changed <- get_parameters("post_disturb")
  parameter_list <- vary_parameter_value(param_changed, variable_interest, vector_value)
  Initial_List <- simulate_predisturb_initial_list()
  
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


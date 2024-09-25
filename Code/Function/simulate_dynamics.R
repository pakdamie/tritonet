#' Force a disturbance (used in conjuction for ODE)
#'
#'This is the function that you need to modify if you want
#'to include a disturbance in deSolve
#'
#' @param num_patch Number of patches that you want to simulate
#' @param frequency How many times are the patches disturbed
#' @param coverage The proportion (0-1) of the patches that is targeted
#' @param species1  How much of the primary vector species survive
#' @param species2 How much of the secondary vector species survive
#' @param end_time total length of the simulation
#' @param randomize is the patches chosen randomized each disturbance?
#'
#' @return a data.frame
#'
#' @examples
force_event_ODE <- function(num_patch, frequency, coverage,
                            species1, species2, end_time, randomize = FALSE) {
  if (randomize == FALSE) {
    # Sample from all the different patches.
    chosen_patch <- matrix(sample(seq(1, num_patch),
      floor(num_patch * coverage),
      replace = FALSE
    ))

    disturbance_df <- data.frame(
      var = c(namer_chosen_compartments(chosen_patch)),
      value = c(
        rep(species1, 2 * length(chosen_patch)),
        rep(species2, 2 * length(chosen_patch))
      ),
      method = c(rep("mult", 4 * length(chosen_patch)))
    )

    event_df <- disturbance_df[rep(
      seq_len(nrow(disturbance_df)),
      length(seq(1, end_time, frequency))
    ), ]

    event_df$time <- rep(seq(1, end_time, frequency),
      each = nrow(disturbance_df)
    )
  } else {
    times_spray <- seq(1, end_time, frequency)

    disturbance_list <- NULL

    for (spray in seq(1, length(times_spray))) {
      chosen_patch_sample <- matrix(sample(seq(1, num_patch),
        floor(num_patch * coverage),
        replace = FALSE
      ))

      disturbance_df <- data.frame(
        var = c(namer_chosen_compartments(chosen_patch_sample)),
        value = c(
          rep(species1, 2 * length(chosen_patch_sample)),
          rep(species2, 2 * length(chosen_patch_sample))
        ),
        method = c(rep("mult", 4 * length(chosen_patch_sample)))
      )

      disturbance_df$time <- times_spray[[spray]]

      disturbance_list[[spray]] <- disturbance_df
    }
    event_df <- do.call(rbind, disturbance_list)
  }
  return(event_df)
}


#' Simulate the primary and secondary vector on the network
#' 
#' This is the main simulation function for the model
#'
#' @param num_patch Number of patches that you want to simulate
#' @param connectance The connectance of the network
#' @param max_distance Default = 20, The total maximum distance
#' @param coverage The proportion of the patches that is targeted
#' @param frequency How many times are the patches disturbed
#' @param species1 How much of the primary vector species survive
#' @param species2 How much of the secondary vector species survive
#' @param initial_values "default" uses the standard setting, but user
#' can input a data.frame.
#' @param parameter_values 'default' uses the standard setting, but
#' user can input a data.frame
#' @param disturbance 'default' - yes, includes disturbance in the model
#' @param disease_on 'default' - no 
#' @param end_length total length of the simulation
#' @param seed the random seed generator 
#' @return
#' @export
#'
#' @examples
#' 
#' 
###Networks of interest

Simulate_model  <- function(user_network,
                               max_distance,
                               coverage,
                               frequency,
                               species1,
                               species2,
                               initial_values = "default",
                               parameter_values = "default",
                               disturbance = 'yes', 
                               disease_on = 'no',
                               simulation_time,
                               intervention_time,
                               seed = NULL,
                               randomize = FALSE){
        
        num_patch <- vcount(user_network)
        
        adjacency_matrix <- Matrix(as_adjacency_matrix(
                user_network,
                type = c("both"),
                attr = "weight",
                names = TRUE,
                sparse = F),
        sparse = TRUE)
        
        if (initial_values == "default"){
                
        initial_y <- c(HS =  rep(1000, num_patch, num_patch),
                        HI = rep(0,num_patch),
                        HR = rep(0,num_patch),
                        PS = sample(seq(1000,2000), num_patch, replace = TRUE),
                        PI = rep(100,num_patch),
                        SS = sample(seq(1000,2000), num_patch, replace = TRUE),
                        SI = rep(100,num_patch))        
                
                
        ###The initial conditions
        #initial_y <- c(HS = rep(1000, num_patch),
        #               HI = rep(5,num_patch),
        #               HR = rep(0,num_patch),
        #               PS = rep(100, num_patch) ,
        #               PI = rep(10,num_patch),
        #               SS = rep(100, num_patch),
        #               SI = rep(5,num_patch))
        }else{
        initial_y <-  initial_values 
        }

        if(parameter_values == "default"){
                
        ###These are daily rates
        parameters_full <- c(
                b_H = 1/27375, #Human birth rate
                b_P = (0.05 + 1e-6), #P.vector birth rate
                b_S = 0.05, #S. vector birth rate
                mu_H = 1/27375, #Human death rate
                mu_P = 0.05, #P. vector death rate
                mu_S = 0.05, #S. vector death rate
                
                a_P = 7, #daily biting rate of the p. vector
                a_S = 7 * (0.55), #daily biting rate of the s.vector
                
                phi_P = 0.132, #transmission probability of p. vector
                phi_S = 0.132 * 0.55, #transmission probability of s. vector
                phi_H  = 0.116, #transmission probability of human
                
                # Recovery rate of the acute phase
                gamma = 1/56,  #recovery rate of infected human
                
                #competition coefficient
                c_PS = 5e-4,#competitition effect of p.vector on s.vector
                c_SP = 1e-6,  #competitition effect of s.vector on p.vector
                
                a_max = 2 * 10^-2 ,
                k = 0.02,
                a_0 = 300
               
                 )
        } else{
                parameters_full <-   parameter_values
        }
        
        if(disease_on == "no"){
                parameters_full["phi_P"] <- 0
                parameters_full["phi_S"] <- 0
                parameters_full["phi_H"] <- 0
        }else{
                print("Note: Disease is spreading!")
        }
        
        
        event_DF <- force_event_ODE(num_patch, frequency, coverage,
                        species1, species2, intervention_time, randomize)
        
        
        # When disturbance is "yes" then ensure that you included
        # the event_DF 
        if(disturbance == "yes"){
                
                results <- data.frame(ode(
                        times = seq(1, simulation_time,1),
                        y =   initial_y  ,
                        func = model_ross_trito_metapopulation,
                        parms = parameters_full,
                        num_patch =  num_patch ,
                        adj_matrix = adjacency_matrix,
                        events = list(data = event_DF),
                        method = 'lsodar'),
                        coverage,
                        frequency)
                
        }
        else {

        results <- data.frame(ode(
                times = seq(1,simulation_time,1),
                y =   initial_y  ,
                func = model_ross_trito_metapopulation,
                parms = parameters_full,
                num_patch =  num_patch ,
                adj_matrix = adjacency_matrix,
                method = 'lsodar'),
                coverage,
                frequency)
                
        }
        return((results))
}

#' Simulate the primary and secondary vector on the network
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
Simulator_function <- function(num_patch,
                               connectance,
                               max_distance = 20,
                               coverage,
                               frequency,
                               species1 = 0.8,
                               species2 = 0.2,
                               initial_values = "default",
                               parameter_values = "default",
                               disturbance = 'yes', 
                               disease_on = 'no',
                               end_length, 
                               seed = 24601){
        
        adjacency_matrix <- simulate_final_adjacency_matrix(seed, num_patch,
                                        connectance,max_distance)
        
        g9 <- graph_from_adjacency_matrix(adjacency_matrix, weighted=TRUE,
                                          mode="plus", diag=FALSE)
        
        initial_vec <- sample(seq(1,100), num_patch, replace = TRUE)
        
        if(initial_values == "default"){
        ###The initial conditions
        initial_y <- c(HS = sample(seq(1,1000),num_patch, replace = TRUE),
                       HI = rep(5,num_patch),
                       HR = rep(0,num_patch),
                       PS =       initial_vec ,
                       PI = rep(10,num_patch),
                       SS =        initial_vec ,
                       SI = rep(5,num_patch))
        }else{
        initial_y <-  initial_values 
        }

        if( parameter_values == "default"){
        parameters_full <- c(
                b_H = 1/27375, #Human birth rate
                b_P = 0.09, #P.vector birth rate
                b_S = 0.09, #S. vector birth rate
                mu_H = 1/27375, #Human death rate
                mu_P =0.09, #P. vector death rate
                mu_S = 0.09, #S. vector death rate
                
                a_P = 4, #biting rate of the p. vector
                a_S = 2, #biting rate of the s.vector
                
                phi_P =  0.0009, #transmission probability of p. vector
                phi_S = 0.00004, #transmission probability of s. vector
                phi_H  = 0.005, #transmission probability of human
                
                # Recovery rate of the acute phae
                gamma = 1/56,  #recovery rate of infected human
                
                #competition coefficient
                c_PS = 1e-4,#competitition effect of p.vector on s.vector
                c_SP = 1e-6,  #competitition effect of s.vector on p.vector
                
                a_max =1,
                k = 1e-2,
                a_0 = 250, 
                
                lambda = 0.01)
        } else{
                parameters_full =    parameter_values
        }
        
        if(disease_on == "no"){
                parameters_full["phi_P"]<- 0
                parameters_full["phi_S"]<-0
                parameters_full["phi_H"]<- 0
        }else{
                print("Disease is spreading!")
        }
        
        event_DF <-force_event_ODE(num_patch, frequency, coverage,
                        species1, species2, end_length)
        
        if(disturbance == "yes"){
                results <- data.frame(ode(
                        times = seq(1,  end_length,1),
                        y =   initial_y  ,
                        func = model_ross_trito_metapopulation,
                        parms = parameters_full,
                        num_patch =  num_patch ,
                        adj_matrix = adjacency_matrix,
                        events = list(data = event_DF),
                        method = 'lsodar'))
                
        }
        else{

        results <- data.frame(ode(
                times = seq(1,end_length,1),
                y =   initial_y  ,
                func = model_ross_trito_metapopulation,
                parms = parameters_full,
                num_patch =  num_patch ,
                adj_matrix = adjacency_matrix,
                method = 'lsodar'))
                
        }
        return(list(results, unique(event_DF$var)))
}

###Test
#results_ode<- Simulator_function(num_patch = 25,connectance = 0.05,
#                 max_distance = 20,
#                  coverage = 0.1, frequency = 10, 
#                  species1 = 0.3, species2 = 0.8,
#                  initial_values = "default",
#                  parameter_values = "default",
#                  disturbance = "yes",         
#                  disease_on = 'no',
#                  end_length = 300)

        

#' Simulate the 2-vector/1 host model on a network
#'
#' This is the main model that we simulate infections on.
#'
#'
#' @param ntime How long should this run for? (numeric)
#' @param param The parameter data.frame (data.frame)
#' @param disturbance_time The time that you disturb the system (numeric)
#' @param mortality_P The proportion of primary vector that survives (numeric)
#' @param mortality_S  The proportion of secondary vector that survives (numeric)
#' @param coverage The proportion of the patches randomly chosen (numeric)
#' @param user_network The network to run the model on (igraph object)
#' @param delta_T The time-step
#' @return A list containing seven lists and R effective
#' @export
#'
#' @examples discrete_trito_model(100, param_df, 50, 0.30, 0.70, 0.5, 
#' network, 1 )
#' 
discrete_trito_model_ALTDISP <- function(ntime, 
                                 param,
                                 disturbance_time,
                                 mortality_P,
                                 mortality_S,
                                 coverage, 
                                 user_network, 
                                 delta_T) {
        
        ### Convert the igraph object as an adjacency matrix 
        
        adjacency_matrix <- as_adjacency_matrix(
                user_network,
                type = "both",
                attr = "weight",
                names = TRUE,
                sparse = FALSE
        )
        
        summed_prob <- rowSums(adjacency_matrix)
        adjacency_matrix_adj <- sweep(adjacency_matrix, MARGIN = 1, summed_prob, `/`)
        
        #How many compartments do we need?
        patch_num <- nrow(adjacency_matrix)
        
        # We create matrices for the different classes that must be tracked
        compartment_label <- c(
                "HS_mat", "HI_mat", "HR_mat",  # Humans (sus/inf/rec)
                "PS_mat", "PI_mat",  # Primary (sus/inf)
                "SS_mat", "SI_mat"   # Secondary(sus/inf)
        )
        
        # Using the above compartment label, we then create an empty matrix where 
        # the row are the time-steps and the columns are the individual patches
        
        for (i in 1:length(compartment_label)) {
                assign(
                        compartment_label[i],
                        matrix(0, nrow = ntime, ncol = patch_num)
                )
        }
        
        #Initial conditions (everyone starts out with the same number)
        HS_mat[1, ] <- rep(1000, patch_num)
        HI_mat[1, ] <- rep(0, patch_num) 
        HR_mat[1, ] <- rep(0, patch_num)
        
        PS_mat[1, ] <- rep(1000, patch_num)
        PI_mat[1, ] <- rep(100, patch_num) 
        SS_mat[1, ] <- rep(1000, patch_num)
        SI_mat[1, ] <- rep(100, patch_num)
        
        #Run the parameter values to load
        b_H <- param["b_H"] # Birth rate of humans
        b_P <- param["b_P"] # Birth rate of primary vectors
        b_S <- param["b_S"] # Birth rate of secondary vectors
        mu_H <- param["mu_H"] # Mortality rate of humans
        mu_V <- param["mu_V"] # Mortality rate of vectors
        c_PS <- param["c_PS"] # Competition coefficient of primary on secondary
        c_SP <- param["c_SP"] # Competition coefficient of secondary on primary
        c_PP <- param["c_PP"] # Competition coefficient of primary on primary
        c_SS <- param["c_SS"] # Competition coefficient of secondary on secondary
        f_P <- param["f_P"]   # Biting rate of primary
        f_S <- param["f_S"]   # Biting rate of secondary
        theta_P <- param["theta_P"] # Transmission probability of primary
        theta_S <- param["theta_S"] # Transmission probability of secondary
        theta_H <- param["theta_H"] # Transmission probability of humans
        gamma <- param["gamma"] # Recovery rate
        d <- param["d"]
        
        for (j in 1:(ntime - 1)) {
                
                ### Retrieve the total populations in patch j
                NH_mat <- HS_mat[j,  ] + HI_mat[j, ] + HR_mat[j, ] # Total humans
                NP_mat <- PS_mat[j, ] + PI_mat[j, ] # Total primary vectors
                NS_mat <- SS_mat[j, ] + SI_mat[j, ] # Total secondary vector population
                
                #The if-else statements prevent situations where we divide by 0 
                # Is the equation finite?
                HS_ratio <- ifelse(!(is.finite(SI_mat[j, ] / NS_mat)),
                                   0, #If not finite
                                   SI_mat[j, ]/ NS_mat #If finite
                )
                
                HP_ratio <- ifelse(
                        !is.finite(PI_mat[j, ] / NP_mat), 
                        0, #If not finite
                        PI_mat[j, ] / NP_mat #If finite
                )
                
                H_ratio <- ifelse(
                        !is.finite(HI_mat[j, ] / NH_mat),
                        0, #If not finite
                        HI_mat[j, ] / NH_mat #If finite
                )
                
                # Primary and secondary will infect humans
                infections_H <- (theta_P * f_P * HP_ratio) + (theta_S * f_S * HS_ratio)  
                
                #Both primary and secondary vectors are infected by humans
                infections_P <- (theta_H * f_P) *  H_ratio # Primary infections
                infections_S <- (theta_H * f_S) *  H_ratio # Secondary infections
                
                #Demographic/infection rates of humans
                HS_Rates <- (b_H * NH_mat) - (infections_H + mu_H) * HS_mat[j, ]
                HI_Rates <- (infections_H * HS_mat[j, ]) - (gamma + mu_H) * HI_mat[j, ]
                HR_Rates <- (gamma * HI_mat[j, ]) - (mu_H * HR_mat[j, ])
                
                #Demographic/infection rates of primary vector
                PS_Rates <- (b_P * NP_mat) - (infections_P + mu_V + 
                                                      (c_SP * NS_mat) + (c_PP * NP_mat)) * PS_mat[j, ] 
                
                PI_Rates <- (infections_P * PS_mat[j, ]) - (mu_V  + 
                                                                    (c_SP * NS_mat) + (c_PP * NP_mat)) * PI_mat[j, ]
                
                #Demographic/infection of secondary vector
                SS_Rates <- (b_S * NS_mat) - (infections_S + mu_V + 
                                                      (c_PS * NP_mat) + (c_SS * NS_mat)) * SS_mat[j, ] 
                
                SI_Rates <- (infections_S * SS_mat[j, ]) - (mu_V +  
                                                                    (c_PS * NP_mat) + (c_SS * NS_mat)) * SI_mat[j, ]
                
                
                # Dispersal of the primary and secondary vectors
                dispersal_PS <- rowSums(sweep(adjacency_matrix, MARGIN =1, (- d * PS_mat[j, ]),`*`)) + 
                        t(adjacency_matrix) %*% (d * PS_mat[j, ]) 
                
                dispersal_PI <-  rowSums(sweep(adjacency_matrix, MARGIN =1, (- d * PI_mat[j, ]),`*`)) + 
                        t(adjacency_matrix) %*% (d * PI_mat[j, ]) 
                
                dispersal_SS <- rowSums(sweep(adjacency_matrix, MARGIN =1, (- d * SS_mat[j, ]),`*`)) + 
                        t(adjacency_matrix) %*% (d * SS_mat[j, ]) 
                
                dispersal_SI <- rowSums(sweep(adjacency_matrix, MARGIN =1, (- d * SI_mat[j, ]),`*`)) + 
                        t(adjacency_matrix) %*% (d * SI_mat[j, ]) 
                
                
                #Let's look what the new population size will be
                total_change_HS <- HS_mat[j, ] + (HS_Rates * delta_T)
                total_change_HI <- HI_mat[j, ] + (HI_Rates* delta_T)
                total_change_HR <- HR_mat[j, ] + (HR_Rates * delta_T)
                total_change_PS <- PS_mat[j, ] + (PS_Rates + dispersal_PS) * delta_T
                total_change_PI <- PI_mat[j, ] + (PI_Rates + dispersal_PI) * delta_T
                total_change_SS <- SS_mat[j, ] + (SS_Rates + dispersal_SS) * delta_T
                total_change_SI <- SI_mat[j, ] + (SI_Rates + dispersal_SI) * delta_T
                
                total_change_HS [total_change_HS  < 0] <- 0
                total_change_HI [total_change_HI  < 0] <- 0
                total_change_HR [total_change_HR  < 0] <- 0
                total_change_PS [total_change_PS  < 0] <- 0
                total_change_PI [total_change_PI  < 0] <- 0
                total_change_SS [total_change_SS  < 0] <- 0
                total_change_SI [total_change_SI  < 0] <- 0
                
                
                #If there is no disturbance, update!
                if (j != disturbance_time) {
                        
                        #If it's negative, make everything 0
                        
                        #We then update these row by row
                        HS_mat[j + 1, ] <- total_change_HS 
                        HI_mat[j + 1, ] <- total_change_HI
                        HR_mat[j + 1, ] <- total_change_HR
                        PS_mat[j + 1, ] <- total_change_PS 
                        PI_mat[j + 1, ] <- total_change_PI 
                        SS_mat[j + 1, ] <- total_change_SS 
                        SI_mat[j + 1, ] <- total_change_SI 
                }
                
                # IF there is a disturbance
                else if (j == disturbance_time) {
                        
                        ### Coverage gives you the column indices of the patches 
                        coverage <- sample(1:patch_num, size = floor(coverage *  patch_num ))
                        
                        HS_mat[j + 1, ] <- total_change_HS 
                        HI_mat[j + 1, ] <- total_change_HI
                        HR_mat[j + 1, ] <- total_change_HR
                        PS_mat[j + 1, ] <- total_change_PS 
                        PI_mat[j + 1, ] <- total_change_PI 
                        SS_mat[j + 1, ] <- total_change_SS 
                        SI_mat[j + 1, ] <- total_change_SI 
                        
                        PS_mat[j+1, coverage] <- (PS_mat[j+1, coverage]) *  mortality_P
                        SS_mat[j+1, coverage] <- (SS_mat[j+1, coverage]) *  mortality_S      
                        
                        PI_mat[j+1, coverage] <- (PI_mat[j+1, coverage]) *  mortality_P
                        SI_mat[j+1, coverage] <- (SI_mat[j+1, coverage]) *  mortality_S     
                }
        }
        
        list_abundance <- list(HS_mat, HI_mat, HR_mat,
                               PS_mat, PI_mat, 
                               SS_mat, SI_mat)
        
        #When the loop has ended- put the matrices into the list
        
        return(list_abundance)
        
}


calculate_R_effective_discrete_patch <- function(parameters, 
                                                 patch_num,
                                                 NH, NP, NS,
                                                 HS, HI, HR,
                                                 PS, PI,
                                                 SS, SI) {
        
        theta_P <- parameters["theta_P"]
        theta_S <- parameters["theta_S"]
        theta_H <- parameters["theta_H"]
        gamma <- parameters["gamma"]
        mu_H <- parameters["mu_H"]
        mu_V <- parameters["mu_V"]
        f_P <- parameters["f_P"]
        f_S <- parameters["f_S"]
        c_SP <- parameters["c_SP"]
        c_PS <- parameters["c_PS"]
        c_SS <- parameters["c_SS"]
        c_PP <- parameters["c_PP"]
        
        human_R0_patch <- NULL
        
        for (k in seq(1, patch_num)){
        
                
          H_P_ratio <- ifelse(is.finite(HS[k]/NP[k]), HS[k]/NP[k], 0)
          H_S_ratio <- ifelse(is.finite(HS[k]/NS[k]), HS[k]/NS[k], 0)
          P_H_ratio <- ifelse(is.finite(PS[k]/NH[k]), PS[k]/NH[k], 0)
          S_H_ratio <- ifelse(is.finite(SS[k]/NH[k]), SS[k]/NH[k], 0)
          
              
          F_mat <- matrix(c(0, (theta_P * f_P * H_P_ratio), (theta_S * f_S * H_S_ratio),
                          (theta_H * f_P * P_H_ratio), 0, 0,
                           (theta_H * f_S *  S_H_ratio), 0, 0), 
                        byrow = TRUE, ncol = 3) 
          

          F_mat[!(is.finite(F_mat))] <- 0
                
          V_mat <- matrix(c((gamma + mu_H), 0, 0,
                            0, mu_V + (c_SP * NS[k]) + (c_PP * NP[k]), 0,
                            0, 0, mu_V + (c_PS * NP[k]) + (c_SS * NS[k])),
                            byrow = TRUE, ncol = 3)
                
          human_R0_patch[[k]] <- max(eigen(F_mat %*% solve(V_mat))$values)
        }
        
        return(do.call(rbind, human_R0_patch))
        
}

#' Simulate the 2-vector/1 host model on a network
#'
#' This is the main model that we simulate infections with.
#'
#'
#' @param ntime How long should this run for?
#' @param param The parameter data.frame 
#' @param disturbance_time The time point that you disturb the system
#' @param mortality_P The proportion of primary vector that survives 
#' @param mortality_S  The proportion of secondary vector that survives 
#' @param coverage The proportion of the patches randomly chosen.
#' @param user_network The network to run the model on
#'
#' @return A list containing lists and R effective
#' @export
#'
#' @examples
#' 
discrete_trito_model <- function(ntime, 
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
 
  #How many compartments do we need 
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
  
  #Initial conditions
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
  d <- param["d"] # Dispersal rate

 
 List_R_effective <- NULL
         
  for (j in 1:(ntime - 1)) {
          
      ### Retrieve the total populations in patch j
      NH_mat <- HS_mat[j,  ] + HI_mat[j, ] + HR_mat[j, ] # Total humans
      NP_mat <- PS_mat[j, ] + PI_mat[j, ] # Total primary vectors
      NS_mat <- SS_mat[j, ] + SI_mat[j, ] # Total secondary vector population

      #The if-else statements prevent situations where we divide by 0 
      HS_ratio <- ifelse(!(is.finite(SI_mat[j, ] / NS_mat)),
        0,
        SI_mat[j, ]/ NS_mat 
      )
      
      HP_ratio <- ifelse(
        !is.finite(PI_mat[j, ] / NP_mat),
        0,
        PI_mat[j, ] / NP_mat
      )
      
      H_ratio <- ifelse(
        !is.finite(HI_mat[j, ] / NH_mat),
        0,
        HI_mat[j, ] / NH_mat
      )
      
      # Primary and secondary will infect humans
      infections_H <- (theta_P * f_P * HP_ratio) + (theta_S * f_S * HS_ratio)  
                      
      #Both primary and secondary vectors are infected by humans
      infections_P <- (theta_H * f_P) *  H_ratio # Primary infections
      infections_S <- (theta_H * f_S) *  H_ratio # Secondary infections

      #Rates of humans
      HS_Rates <- (b_H * NH_mat) - (infections_H + mu_H) * HS_mat[j, ]
      HI_Rates <- (infections_H * HS_mat[j, ]) - (gamma + mu_H) * HI_mat[j, ]
      HR_Rates <- (gamma * HI_mat[j, ]) - (mu_H * HR_mat[j, ])

      #Demographic/infection of primary vector
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
      dispersal_PS <- 
         -(t(adjacency_matrix) %*% (d * PS_mat[j, ]) / rowSums(adjacency_matrix))+
        (adjacency_matrix %*% (d * PS_mat[j, ]) / rowSums(adjacency_matrix))
      
      dispersal_PI <- -(t(adjacency_matrix) %*% (d * PI_mat[j, ]) / rowSums(adjacency_matrix))+
        (adjacency_matrix %*% (d * PI_mat[j, ]) / rowSums(adjacency_matrix))
      
      dispersal_SS <- -(t(adjacency_matrix) %*% (d * SS_mat[j, ]) / rowSums(adjacency_matrix))+
        (adjacency_matrix %*% (d * SS_mat[j, ]) / rowSums(adjacency_matrix))
      
      dispersal_SI <- 
              -(t(adjacency_matrix) %*% (d * SI_mat[j, ]) / rowSums(adjacency_matrix))+
        (adjacency_matrix %*% (d * SI_mat[j, ]) / rowSums(adjacency_matrix))
      
      ###Calculate the R0 before thing changes.
      List_R_effective [[j]] <- calculate_R_effective_discrete_patch(
        param, 
        patch_num,
        NH = NH_mat, 
        NP = NP_mat, 
        NS = NS_mat,
        HS = HS_mat[j,],
        HI = HI_mat[j ],
        PS = PS_mat[j, ],
        PI = PI_mat[j, ],
        SS = SS_mat[j, ],
        SI = SI_mat[j, ])
      
      
      #Let's look what the new population size will be
      total_change_HS <- HS_mat[j, ] + HS_Rates
      total_change_HI <- HI_mat[j, ] + HI_Rates
      total_change_HR <- HR_mat [j, ] + HR_Rates
      
      total_change_PS <- PS_mat[j, ] + (PS_Rates + dispersal_PS)*delta_T
      total_change_PI <- PI_mat[j, ] + (PI_Rates + dispersal_PI)*delta_T
      total_change_SS <- SS_mat[j, ] + (SS_Rates + dispersal_SS)*delta_T
      total_change_SI <- SI_mat[j, ] + (SI_Rates + dispersal_SI)*delta_T

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
  
  return(list(list_abundance,
              List_R_effective))
  
  
}


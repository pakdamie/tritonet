
#' Discrete-time metapopulation model with human hosts and two vectors
#'
#'
#'
#' @param ntime How long to simulate the model for
#' @param param The data.frame to feed into the model 
#' @param adjacency_matrix The adjacency matrix of the network
#' 
#' @return A list with 7 elements: HS, HI, HR, PS, PI, SS, SI
#' @export
#'
#' @examples discrete_trito_model(100, param_df, adjacency_matrix)

discrete_trito_model <- function(ntime,
                                 param,
                                 adjacency_matrix, 
                                 disturbance_time,
                                 mortality_P,
                                 mortality_S,
                                 coverage
                                 ) {
   
  # The set up of the model
  # The number of patches in the adjacency matrix should be the nrow/ncol        
  patch_num <- nrow(adjacency_matrix)
  
  # We need to create matrices for the different classes that must 
  # be tracked
  compartment_label <- c(
    "HS_mat", "HI_mat", "HR_mat",
    "PS_mat", "PI_mat",
    "SS_mat", "SI_mat"
  )

  # Using the above compartment label, we then create an 
  # empty matrix where the row are the time-steps and
  # the columns are the individual patches
  
  for (i in 1:length(compartment_label)) {
    assign(
      compartment_label[i],
      matrix(0, nrow = ntime, ncol = patch_num)
    )
  }
  
  HS_mat[1, ] <- rep(1000, patch_num )
  HI_mat[1, ] <- rep(50, patch_num )
  HR_mat[1, ] <- rep(0, patch_num)
  
  PS_mat[1, ] <- sample(300:500, patch_num) 
  PI_mat[1, ] <- sample(0:10, patch_num) 
  SS_mat[1, ] <- sample(300:500, patch_num)
  SI_mat[1, ] <- sample(0:10, patch_num)
  
  #The parameters
  b_H <- param["b_H"] #birth rate of humans
  b_P <- param["b_P"] #birth rate of primary vectors
  b_S <- param["b_S"] #birth rate of secondary vectors
  mu_H <- param["mu_H"] #mortality rate of human
  mu_V <- param["mu_V"] #mortality rate of vectors
  f_P <- param["f_P"] #biting rate of primary
  f_S <- param["f_S"] #biting rate of secondary
  theta_P <- param["theta_P"] #transmission probability of primary
  theta_S <- param["theta_S"] #transmission probability of secondary
  theta_H <- param["theta_H"] #transmission probability of humans
  gamma <- param["gamma"] #recovery rate
  c_PS <- param["c_PS"] #competition coefficient of primary on secondary
  c_SP <- param["c_SP"] #competition coefficient of secondary on primary
  d <- param["d"] #dispersal rate
  
  for (j in 1 : (ntime - 1)) {
          
      ### Retrieve the total population in patch to calculate frequency-dependent
      ### infection.
          
      NH_mat <- HS_mat[j,  ] + HI_mat[j, ] + HR_mat[j, ] # Total humans
      NP_mat <- PS_mat[j, ] + PI_mat[j, ] # Total primary vectors
      NS_mat <- SS_mat[j, ] + SI_mat[j, ] # Total secondary vector population

      # For each class, calculate the dispersal in and out of the patches. The
      # first term is the dispersal out AND second term is the incoming vectors
      # from neighboring patches and normalize them. 
      
      dispersal_PS <- c(-(d * PS_mat[j, ]) + (adjacency_matrix %*% (d * PS_mat[j, ]) /
        rowSums(adjacency_matrix)))

      dispersal_PI <- c(-(d * PS_mat[j, ]) + (adjacency_matrix %*% (d * PI_mat[j, ]) /
        rowSums(adjacency_matrix)))
      
      dispersal_SS <- c(-(d * PS_mat[j, ]) + (adjacency_matrix %*% (d * SS_mat[j, ]) /
        rowSums(adjacency_matrix)))
      
      dispersal_SI <- c(-(d * PS_mat[j, ]) + (adjacency_matrix %*% (d * SI_mat[j, ]) /
        rowSums(adjacency_matrix)))
      
      # Susceptible humans in patches are infected by either primary
      # or secondary vector (assuming simple mass action)
     
      #The if else statements prevent situations where the force
      # of infection is divided by 0 (for example, if the total vector
      # population is wiped out in one patch...)
      HS_ratio <- ifelse(!(is.finite((SI_mat[j,]/ NS_mat))),
                         0,SI_mat[j,]/NS_mat)
      
      HP_ratio <- ifelse(!(is.finite((PI_mat[j,]/ NP_mat))),
                         0,PI_mat[j,]/NP_mat)
      
      
      H_ratio <-ifelse(!(is.finite((HI_mat[j,]/ NH_mat))),
                       0,HI_mat[j,]/NH_mat)
      
      
      infections_H <- (theta_P * f_P) * HP_ratio 
                      (theta_S * f_S) * HS_ratio  
                      
      
      # Both primary and secondary vectors are infected by humans
      infections_P <- (theta_H * f_P) *  H_ratio # Primary infections
      infections_S <- (theta_H * f_S) *  H_ratio# Secondary infections

      # Rates of transmission for humans
      HS_Rates <- (b_H * NH_mat) - (infections_H+ mu_H) * c(HS_mat[j, ])
      HI_Rates <- (infections_H * HS_mat[j, ]) - (gamma + mu_H) * HI_mat[j, ]
      HR_Rates <- gamma * (HI_mat[j, ]) - (mu_H * HR_mat[j, ])

      # Rates of transmission for primary vectors
      PS_Rates <- (b_P * NP_mat) + dispersal_PS - (infections_P + mu_V + (c_SP * NS_mat)) * PS_mat[j, ] 
      PI_Rates <- dispersal_PI + (infections_P * PS_mat[j, ]) - (mu_V  + (c_SP * NS_mat)) * PI_mat[j, ]

      # ... and Secondary vectors
      SS_Rates <- (b_S * NS_mat) + dispersal_SS - (infections_S + mu_V * (c_PS * NP_mat)) * SS_mat[j, ] 
      SI_Rates <- dispersal_SI + (infections_S * SS_mat[j, ]) - (mu_V +  (c_PS * NP_mat)) * SI_mat[j, ]

    
      # Let's look what the new population size will be
      total_change_HS <- HS_mat[j, ] + HS_Rates
      total_change_HI <- HI_mat[j, ] + HI_Rates
      total_change_HR <- HR_mat [j, ] + HR_Rates
      total_change_PS <- PS_mat[j, ] + PS_Rates
      total_change_PI <- PI_mat[j, ] + PI_Rates
      total_change_SS <- SS_mat[j, ] + SS_Rates
      total_change_SI <- SI_mat[j, ] + SI_Rates
      
      #If it's negative, make everything 0
      total_change_HS [total_change_HS  < 0] <- 0
      total_change_HI [total_change_HI  < 0] <- 0
      total_change_HR [total_change_HR  < 0] <- 0
      total_change_PS [total_change_PS  < 0] <- 0
      total_change_PI [total_change_PI  < 0] <- 0
      total_change_SS [total_change_SS  < 0] <- 0
      total_change_SI [total_change_SI  < 0] <- 0

      if(j != disturbance_time){
      #We then update these row by row
      HS_mat[j + 1, ] <- total_change_HS
      HI_mat[j + 1, ] <- total_change_HI
      HR_mat[j + 1, ] <- total_change_HR
      PS_mat[j + 1, ] <- total_change_PS
      PI_mat[j + 1, ] <- total_change_PI
      SS_mat[j + 1, ] <- total_change_SS
      SI_mat[j + 1, ] <- total_change_SI
      }
      
      else if (j == disturbance_time){
      
      coverage <- sample(1: patch_num , size = floor(coverage *  patch_num ))
                        
      total_change_PS [coverage] <- PS_mat[j,coverage ] *  mortality_P
      total_change_PI [coverage] <- PI_mat[j, coverage ] *  mortality_P
      total_change_SS[coverage] <- SS_mat[j, coverage] * mortality_S
      total_change_SI [coverage] <- SI_mat[j, coverage] * mortality_S
     
      
      #We then update these row by row
      HS_mat[j + 1, ] <- total_change_HS
      HI_mat[j + 1, ] <- total_change_HI
      HR_mat[j + 1, ] <- total_change_HR
      PS_mat[j + 1, ] <- total_change_PS
      PI_mat[j + 1, ] <- total_change_PI
      SS_mat[j + 1, ] <- total_change_SS
      SI_mat[j + 1, ] <- total_change_SI
              
             
      }
  }

  
  
    
  list_abundance <- list(HS_mat, HI_mat, HR_mat,
                         PS_mat, PI_mat, 
                         SS_mat, SI_mat)
  


  #When the loop has ended- put the matrices into the list
  return(list(
    HS_mat, HI_mat, HR_mat,
    PS_mat, PI_mat, 
    SS_mat, SI_mat))
  
  
}


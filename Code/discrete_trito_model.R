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
        
        
        humanR0_patch <- NULL
        
        for (k in seq(1, patch_num)){
        
                
        H_P_ratio <- ifelse(is.finite(HS[k]/NP[k]), HS[k]/NP[k], 0)
        H_S_ratio <- ifelse(is.finite(HS[k]/NS[k]), HS[k]/NS[k], 0)
        P_H_ratio <- ifelse(is.finite(PS[k]/NH[k]), PS[k]/NH[k], 0)
        S_H_ratio <- ifelse(is.finite(SS[k]/NH[k]), SS[k]/NH[k], 0)
          
              
                
                F_mat <- matrix(c(0, (theta_P * f_P * H_P_ratio), (theta_S * f_S * H_S_ratio),
                                  theta_H * f_P * P_H_ratio, 0, 0,
                                  theta_H * f_S *  S_H_ratio, 0, 0), 
                                  byrow = TRUE, ncol = 3) 
        

                F_mat[!(is.finite(F_mat))] <- 0
                
                V_mat <- matrix(c(gamma + mu_H, 0, 0,
                                  0, mu_V + (c_SP * NS[k]), 0,
                                  0, 0, mu_V + (c_PS * NP[k])),
                                  ncol = 3)
                
                humanR0_patch[[k]] <- max(eigen(F_mat %*% solve(V_mat))$values)
        }
        
        
      
        return(do.call(rbind, humanR0_patch))
        
}

discrete_trito_model <- function(ntime, 
                                 param,
                                 disturbance_time,
                                 mortality_P,
                                 mortality_S,
                                 coverage, 
                                 user_network) {
   
 ### We convert the network to an adjacency matrix 
 adjacency_matrix <- as_adjacency_matrix(
                                user_network,
                                type = c("both"),
                                attr = "weight",
                                names = TRUE,
                                sparse = FALSE)    
 
  #The number of patches in the adjacency matrix should be the nrow/ncol 
  patch_num <- nrow(adjacency_matrix)
  
  # We need to create matrices for the different classes that must 
  # be tracked
  compartment_label <- c("HS_mat", "HI_mat", "HR_mat",
                         "PS_mat", "PI_mat",
                         "SS_mat", "SI_mat")

  # Using the above compartment label, we then create an empty matrix where 
  #the row are the time-steps and the columns are the individual patches
  
  for (i in 1:length(compartment_label)) {
    assign(compartment_label[i],
      matrix(0, nrow = ntime, 
             ncol = patch_num)
    )
  }
  
  #We set the initial values of the individuals in the patches
  HS_mat[1, ] <- rep(2000, patch_num)
  HI_mat[1, ] <- rep(0, patch_num)
  HR_mat[1, ] <- rep(0, patch_num)
  
  PS_mat[1, ] <- rep(250, patch_num) 
  PI_mat[1, ] <- rep(10, patch_num) 
  SS_mat[1, ] <- rep(250, patch_num)
  SI_mat[1, ] <- rep(10, patch_num)
  
  #Run the parameter values to load
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
  k_p <- param["k_p"]
  k_s <- param["k_s"]
 
 List_R_effective <- NULL
         
  for (j in 1:(ntime - 1)) {
          
      ### Retrieve the total population in patch to calculate frequency-dependent
      ### infection.
      NH_mat <- HS_mat[j,  ] + HI_mat[j, ] + HR_mat[j, ] # Total humans
      NP_mat <- PS_mat[j, ] + PI_mat[j, ] # Total primary vectors
      NS_mat <- SS_mat[j, ] + SI_mat[j, ] # Total secondary vector population

      #The if else statements prevent situations where the force of
      #infection
      HS_ratio <- ifelse(!(is.finite((SI_mat[j,]/ NS_mat))),
                         0, SI_mat[j, ] / NS_mat)
      
      HP_ratio <- ifelse(!(is.finite((PI_mat[j,] / NP_mat))),
                         0, PI_mat[j, ] / NP_mat)
      
      H_ratio <-ifelse(!(is.finite((HI_mat[j,] / NH_mat))),
                       0, HI_mat[j, ] / NH_mat)
      
      # Primary and secondary will infect humans
      infections_H <- (theta_P * f_P * HP_ratio) + (theta_S * f_S * HS_ratio)  
                      
      #Both primary and secondary vectors are infected by humans
      infections_P <- (theta_H * f_P) *  H_ratio # Primary infections
      infections_S <- (theta_H * f_S) *  H_ratio # Secondary infections

      #Rates of humans
      HS_Rates <- (b_H * NH_mat) - (infections_H + mu_H) * HS_mat[j, ]
      HI_Rates <- (infections_H * HS_mat[j, ]) - (gamma + mu_H) * HI_mat[j, ]
      HR_Rates <- (gamma * HI_mat[j, ]) - (mu_H * HR_mat[j, ])

      #Rates of primary vector
      PS_Rates <- (b_P * NP_mat * (1 - NP_mat/k_p)) - (infections_P + mu_V + (c_SP * NS_mat)) * PS_mat[j, ] 
      PI_Rates <- (infections_P * PS_mat[j, ]) - (mu_V  + (c_SP * NS_mat)) * PI_mat[j, ]

      #Rates of secondary vector
      SS_Rates <- (b_S * NS_mat * (1 - NS_mat/k_s)) - (infections_S + mu_V + (c_PS * NP_mat)) * SS_mat[j, ] 
      SI_Rates <- (infections_S * SS_mat[j, ]) - (mu_V +  (c_PS * NP_mat)) * SI_mat[j, ]
      
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
             
     
      dispersal_PS <- -(d * PS_mat[j, ]) + (adjacency_matrix %*% (d * PS_mat[j, ]) /
                                                    rowSums(adjacency_matrix))
      
      dispersal_PI <- -(d * PI_mat[j, ]) + (adjacency_matrix %*% (d * PI_mat[j, ]) /
                                                    rowSums(adjacency_matrix))
      
      dispersal_SS <- -(d * SS_mat[j, ]) + (adjacency_matrix %*% (d * SS_mat[j, ]) /
                                                    rowSums(adjacency_matrix))
      
      dispersal_SI <- -(d * SI_mat[j, ]) + (adjacency_matrix %*% (d * SI_mat[j, ]) /
                                                    rowSums(adjacency_matrix))
      
      #Let's look what the new population size will be
      total_change_HS <- HS_mat[j, ] + HS_Rates
      total_change_HI <- HI_mat[j, ] + HI_Rates
      total_change_HR <- HR_mat [j, ] + HR_Rates
      
      total_change_PS <- PS_mat[j, ] + PS_Rates + dispersal_PS
      total_change_PI <- PI_mat[j, ] + PI_Rates + dispersal_PI
      total_change_SS <- SS_mat[j, ] + SS_Rates + dispersal_SS
      total_change_SI <- SI_mat[j, ] + SI_Rates + dispersal_SI
      
      total_change_HS [total_change_HS  < 0] <- 0
      total_change_HI [total_change_HI  < 0] <- 0
      total_change_HR [total_change_HR  < 0] <- 0
      total_change_PS [total_change_PS  < 0] <- 0
      total_change_PI [total_change_PI  < 0] <- 0
      total_change_SS [total_change_SS  < 0] <- 0
      total_change_SI [total_change_SI  < 0] <- 0
      

      #If there is no disturbance
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
      
     
      # For each class, calculate the dispersal in and out of the patches. The
      # first term is the dispersal out AND second term is the incoming vectors
      # from neighboring patches and normalize them. 
      
      # Susceptible humans in patches are infected by either primary
      # or secondary vector (assuming simple mass action)

      }
      
      else if (j == disturbance_time) {
      
      coverage <- sample(1: patch_num , size = floor(coverage *  patch_num ))
                        
      print(coverage)
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
  return(list(list_abundance,List_R_effective))
  
  
}


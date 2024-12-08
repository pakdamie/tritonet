calculate_R_effective_discrete_patch <- function(
  parameters, lists, disturbance_time) {
        
  theta_P <- parameters["theta_P"]
  theta_M <- parameters["theta_M"]
  theta_H <- parameters["theta_H"]
  gamma <- parameters["gamma"]
  mu_H <- parameters["mu_H"]
  f_P <- parameters["f_P"]
  f_M <- parameters["f_M"]
  c_MP <- parameters["c_MP"]
  c_PM <- parameters["c_PM"]
  c_MM <- parameters["c_MM"]
  c_PP <- parameters["c_PP"]
        
  HS <- lists[[1]]
  HI <- lists[[2]]
  HR <- lists[[3]]
  
  PS <- lists[[4]]
  PI <- lists[[5]]
  
  MS <- lists[[6]]
  MI <- lists[[7]]
  
  patch_num <- ncol(HS)
  
  NH <- HS + HI + HR #Total human population
  NP <- PS + PI #Total primary population
  NM <- MS + MI #Total secondary population

  ###We're more interested in the RE slightly before disturbance,
  ###during disturbance, and recovery after disturbance
  interest <- seq(disturbance_time - 100, disturbance_time + 400)
  
  full_time <- NULL #(All patches)
 
  for (time in seq_along(interest)){
  
  interest_time <- interest[time]  
  patch_R0_time <- NULL #For a specific patch
  
  for (k in 1:(patch_num)){
                
  H_P_ratio <- HS[interest_time,k]/NH[interest_time,k]
  H_M_ratio <- HS[interest_time,k]/NH[interest_time,k]
  P_H_ratio <- PS[interest_time,k]/NH[interest_time,k]
  M_H_ratio <- MS[interest_time,k]/NH[interest_time,k]

  F_mat <- Matrix(
    c(0, (theta_P * f_P * H_P_ratio), (theta_M * f_M * H_M_ratio),
    (theta_H * f_P * P_H_ratio), 0, 0,
    (theta_H * f_M *  M_H_ratio), 0, 0), 
    byrow = TRUE, ncol = 3,sparse = TRUE) 
                
  H_out <- 1/(gamma + mu_H)
  P_Vout <- (c_MP * NM[interest_time,k]) + (c_PP * NP[interest_time,k])
  M_Vout <-(c_PM * NP[interest_time,k]) + (c_MM * NM[interest_time,k])
  
  P_Vout_Inv <- ifelse(P_Vout == 0, 0, 1/P_Vout)
  M_Vout_Inv <- ifelse(M_Vout == 0, 0, 1/M_Vout)
  
  ### An inverse of a diagonal matrix is the inverse of its diagonal elements
  V_mat <- Matrix(
            c(H_out, 0, 0,
            0, P_Vout_Inv, 0,
            0, 0,  M_Vout_Inv),
            byrow = TRUE, 
            ncol = 3,
            sparse = TRUE)

    
    FV <- (F_mat %*%  V_mat)
    
    ##The RE is the maximum eigenvalue
    RE = max(eigen(FV)$values)
    primary_contrib = sum(FV[,2])
    secondary_contrib = sum(FV[,3])
    
    
    patch_R0_time[[k]] <-  cbind(RE = RE,
                                patch_num = k,
                                time = interest_time,
                                primary_contrib = primary_contrib ,
                                secondary_contrib = secondary_contrib)
  }
  
  full_time[[time]] <- do.call(rbind, patch_R0_time)
  

}


full_time_f <- do.call(rbind.data.frame, full_time)




return(full_time_f)

}


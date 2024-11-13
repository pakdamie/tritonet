calculate_R_effective_discrete_patch <- function(
  parameters, lists, disturbance_time) {
        
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
        
  HS <- lists[[1]]
  HI <- lists[[2]]
  HR <- lists[[3]]
  
  PS <- lists[[4]]
  PI <- lists[[5]]
  
  SS <- lists[[6]]
  SI <- lists[[7]]
  
  
  patch_num <- ncol(HS)
  
  NH <- HS + HI + HR
  NP <- PS + PI
  NS <- SS + SI
  

  interest <- seq(disturbance_time - 4, disturbance_time + 4)
  
full_time <- NULL  
for (time in seq_along(interest)){
  
  interest_time <- interest[time]  
  patch_R0_time <- NULL
  
  for (k in 1:(patch_num)){
                
  H_P_ratio <- ifelse(is.finite(HS[interest_time,k]/NH[interest_time,k]), HS[interest_time,k]/NH[interest_time,k], 0)
  H_S_ratio <- ifelse(is.finite(HS[interest_time,k]/NH[interest_time,k]), HS[interest_time,k]/NH[interest_time,k], 0)
  
  P_H_ratio <- ifelse(is.finite(PS[interest_time,k]/NH[interest_time,k]), PS[interest_time,k]/NH[interest_time,k], 0)
  S_H_ratio <- ifelse(is.finite(SS[interest_time,k]/NH[interest_time,k]), SS[interest_time,k]/NH[interest_time,k], 0)
                
                
  F_mat <- Matrix(c(
    0, (theta_P * f_P * H_P_ratio), (theta_S * f_S * H_S_ratio),
    (theta_H * f_P * P_H_ratio), 0, 0,
    (theta_H * f_S *  S_H_ratio), 0, 0), 
    byrow = TRUE, ncol = 3,sparse = TRUE) 
                
                
    F_mat[!(is.finite(F_mat))] <- 0
                
    V_mat <- Matrix(c((gamma + mu_H), 0, 0,
            0, mu_V + (c_SP * NS[interest_time,k]) + (c_PP * NP[interest_time,k]), 0,
            0, 0, mu_V + (c_PS * NP[interest_time,k]) + (c_SS * NS[interest_time,k])),
            byrow = TRUE, ncol = 3,sparse = TRUE)

           
   if(det(V_mat) == 0){
           RE = 0
           primary_contrib = 0
           secondary_contrib = 0
   }else{
            
    FV <- (F_mat %*%  solve(V_mat))
    RE = max(eigen(FV)$values)
    primary_contrib = sum(FV[1,2])
    secondary_contrib = sum(FV[1,3])
    }
    
    patch_R0_time[[k]] <-  cbind(RE = RE,
                                patch_num = k,
                                time = interest_time,
                                primary_contrib = primary_contrib ,
                                secondary_contrib = secondary_contrib)
        }
  full_time[[time]] <- do.call(rbind, patch_R0_time)
  

}


full_time_f <- do.call(rbind.data.frame, full_time)
maxRE <- full_time_f[which.max(full_time_f$RE),]

max_CV <- max(by(full_time_f, full_time_f$time, function(x)  CV_RE = sd(x$RE)/mean(x$RE)))
maxRE$CV <- max_CV


return(list(full_time_f ,maxRE))

}


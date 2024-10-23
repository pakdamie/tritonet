calculate_R_effective_discrete_patch <- function(
  parameters, patch_num, lists, disturbance_time) {
        
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
  
  NH <- HS + HI + HR
  NP <- PS + PI
  NS <- SS + SI
  

  interest <- seq(disturbance_time -5, disturbance_time + 5)
  
  
full_time <- NULL  
for (t in seq_along(interest)){
  
  interest_time <- interest[t]  
  
  patch_R0_time <- NULL
  
  for (k in 1:(patch_num)){
                
  H_P_ratio <- ifelse(is.finite(HS[interest_time,k]/NP[interest_time,k]), HS[interest_time,k]/NP[interest_time,k], 0)
  H_S_ratio <- ifelse(is.finite(HS[interest_time,k]/NS[interest_time,k]), HS[interest_time,k]/NS[interest_time,k], 0)
  
  
  
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

           
    FV <- F_mat %*%  solve(V_mat)
    
    patch_R0_time[[k]] <- cbind(RE = max(eigen(FV)$values),
                                HI_cont = sum(FV[,1]),
                                PI_cont = sum(FV[,2]),
                                SI_cont = sum(FV[,3]),
                                patch_num = k,
                                time = interest_time)
        }
  full_time[[t]] <- do.call(rbind, patch_R0_time)
  

}





full_time_f <- do.call(rbind.data.frame, full_time)
maxRE <- full_time_f[which.max(full_time_f$RE),]

max_CV <- max(by(full_time_f, full_time_f$time, function(x)  CV_RE = sd(x$RE)/mean(x$RE)))
maxRE$CV <- max_CV


return(list(full_time_f,maxRE))

}


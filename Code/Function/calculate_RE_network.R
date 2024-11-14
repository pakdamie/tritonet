calculate_R_effective_discrete_net <- function(parameters, lists, disturbance_time) {
        
        # Extract parameters
        theta_P <- parameters["theta_P"]
        theta_M <- parameters["theta_M"]
        theta_H <- parameters["theta_H"]
        gamma <- parameters["gamma"]
        mu_H <- parameters["mu_H"]
        mu_V <- parameters["mu_V"]
        f_P <- parameters["f_P"]
        f_M <- parameters["f_M"]
        c_MP <- parameters["c_MP"]
        c_PM <- parameters["c_PM"]
        c_MM <- parameters["c_MM"]
        c_PP <- parameters["c_PP"]
        
        # Extract lists
        HS <- lists[[1]]
        HI <- lists[[2]]
        HR <- lists[[3]]
        
        PS <- lists[[4]]
        PI <- lists[[5]]
        
        MS <- lists[[6]]
        MI <- lists[[7]]
        
        patch_num <- ncol(HS)
        
        # Population counts
        NH <- HS + HI + HR
        NP <- PS + PI
        NS <- MS + MI
        
        # Time window of interest
        interest <- seq(disturbance_time - 5, disturbance_time + 50)
        
        full_time <- list()
        
        for (t in seq_along(interest)) {
                
          interest_time <- interest[t]
                
          # Ratios and infection forces
          
          H_P_ratio <- HS[interest_time,]/NH[interest_time,]
          H_M_ratio <- HS[interest_time,]/NH[interest_time,]
          P_H_ratio <- PS[interest_time,]/NH[interest_time,]
          M_H_ratio <- MS[interest_time,]/NH[interest_time,]
          
                
          P_H_INF <- theta_P * f_P * H_P_ratio
          M_H_INF <- theta_M * f_M * H_M_ratio
          H_P_INF <- theta_H * f_P * P_H_ratio
          H_M_INF <- theta_H * f_M * M_H_ratio
                
          # Initialize matrices
          F_mat <- matrix(0, nrow = length(P_H_INF) * 3, ncol = length( P_H_INF ) * 3, byrow = TRUE)
          V_mat <- matrix(0, nrow = length(P_H_INF) * 3, ncol = length( P_H_INF ) * 3, byrow = TRUE)
                
          for (k in seq_along(P_H_INF)) {
            index <- (3 * (k - 1) + 1) : (3 * k)
                        
           F_mat[index, index] <- matrix(c(
                                0, P_H_INF[k], M_H_INF[k],
                                H_P_INF[k], 0, 0,
                                H_M_INF[k], 0, 0
                        ), ncol = 3, byrow = TRUE)
        }
                

        for (l in seq_along(P_H_INF)) {
          index <- (3 * (l - 1) + 1):(3 * l)
                        
            H_LEAVE <- 1/gamma + mu_H
            P_Vout <- (c_MP * NM[interest_time,]) + (c_PP * NP[interest_time,])
            M_Vout<-  (c_PM * NP[interest_time,]) + (c_MM * NM[interest_time,])
        
            P_Vout_Inv <- ifelse(P_Vout == 0, 0, 1/P_Vout)
            M_Vout_Inv <- ifelse(M_Vout == 0, 0, 1/M_Vout)
                    
                 V_mat[index, index] <- 
                  matrix(c(H_LEAVE, 0, 0,
                           0,  P_Vout_Inv[l], 0,
                           0, 0, M_Vout_Inv[l]
                        ), ncol = 3, byrow = TRUE)
                }
                
                FV <- F_mat %*% (V_mat)
                
          full_time[[t]] <- cbind(RE = max(Re(eigen(FV)$values)),time = interest_time)
        }
                
        
        full_time_f <- do.call(rbind.data.frame, full_time)

        return(full_time_f)
}

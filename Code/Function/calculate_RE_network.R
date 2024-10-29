calculate_R_effective_discrete_net <- function(parameters, lists, disturbance_time) {
        
        # Extract parameters
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
        
        # Extract lists
        HS <- lists[[1]]
        HI <- lists[[2]]
        HR <- lists[[3]]
        PS <- lists[[4]]
        PI <- lists[[5]]
        SS <- lists[[6]]
        SI <- lists[[7]]
        
        patch_num <- ncol(HS)
        
        # Population counts
        NH <- HS + HI + HR
        NP <- PS + PI
        NS <- SS + SI
        
        # Time window of interest
        interest <- seq(disturbance_time - 5, disturbance_time + 5)
        
        full_time <- list()
        
        for (t in seq_along(interest)) {
                
          interest_time <- interest[t]
                
          # Ratios and infection forces
          H_P_ratio <- ifelse(is.finite(HS[interest_time,] / NP[interest_time,]), HS[interest_time,] / NP[interest_time,], 0)
          H_S_ratio <- ifelse(is.finite(HS[interest_time,] / NS[interest_time,]), HS[interest_time,] / NS[interest_time,], 0)
          P_H_ratio <- ifelse(is.finite(PS[interest_time,] / NH[interest_time,]), PS[interest_time,] / NH[interest_time,], 0)
          S_H_ratio <- ifelse(is.finite(SS[interest_time,] / NH[interest_time,]), SS[interest_time,] / NH[interest_time,], 0)
                
          P_H_INF <- theta_P * f_P * H_P_ratio
          S_H_INF <- theta_S * f_S * H_S_ratio
          H_P_INF <- theta_H * f_P * P_H_ratio
          H_S_INF <- theta_H * f_S * S_H_ratio
                
          # Initialize matrices
          F_mat <- matrix(0, nrow = length(H_S_INF) * 3, ncol = length(H_S_INF) * 3, byrow = TRUE)
          V_mat <- matrix(0, nrow = length(H_S_INF) * 3, ncol = length(H_S_INF) * 3, byrow = TRUE)
                
          for (k in seq_along(P_H_INF)) {
            index <- (3 * (k - 1) + 1) : (3 * k)
                        
           F_mat[index, index] <- matrix(c(
                                0, P_H_INF[k], S_H_INF[k],
                                H_P_INF[k], 0, 0,
                                H_S_INF[k], 0, 0
                        ), ncol = 3, byrow = TRUE)
                }
                
                F_mat[!is.finite(F_mat)] <- 0
                
                for (l in seq_along(P_H_INF)) {
                  index <- (3 * (l - 1) + 1):(3 * l)
                        
                    H_LEAVE <- gamma + mu_H
                    P_LEAVE <- mu_V + (c_SP * NS[interest_time,]) + (c_PP * NP[interest_time,])
                    S_LEAVE <- mu_V + (c_PS * NP[interest_time,]) + (c_SS * NS[interest_time,])
                        
                        V_mat[index, index] <- matrix(c(
                                H_LEAVE, 0, 0,
                                0, P_LEAVE[l], 0,
                                0, 0, S_LEAVE[l]
                        ), ncol = 3, byrow = TRUE)
                }
                
                FV <- F_mat %*% solve(V_mat)
                
                full_time[[t]] <- cbind(RE = max(eigen(FV)$values),time = interest_time)
        }
                
        
        full_time_f <- do.call(rbind.data.frame, full_time)
        maxRE <- full_time_f[which.max(full_time_f$RE),]

        max_CV <-  sd(full_time_f $RE) / mean(full_time_f $RE)
        maxRE$CV <- max_CV
        
        return(list(maxRE,  full_time_f ))
}

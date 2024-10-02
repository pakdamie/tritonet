calculate_R_effective_discrete_patch <- function(discrete_lists, 
                                                 parameters, 
                                                 adjacency_matrix,
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
        

        F_mat <- matrix(c(0,theta_P * f_P * HS/NP, theta_S * f_S * HS/NS,
                         theta_H * f_P * PS/NH, 0, 0,
                        theta_H * f_S * SS/NH, 0, 0), byrow = TRUE,ncol =3) 
        
        V_mat <- matrix(c(gamma + mu_H, 0,0,
                               0, mu_V, 0,
                               0, 0, mu_V),ncol =3)
        
  
       max(eigen(F_mat %*% solve(V_mat))$values)
        
                
        
        
        }
        

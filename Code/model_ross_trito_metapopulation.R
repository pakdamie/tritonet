###Modeling metapopulation 

model_ross_trito_metapopulation <- function(t, 
                                            state, 
                                            param, 
                                            num_patch, 
                                            adj_matrix) {
        
        with(as.list(c(state, param)), {

        ###The metapopulation
        H_S = matrix(state[1:num_patch], ncol = 1)
        H_I = matrix(state[(num_patch+1):(2*num_patch)], ncol = 1)
        H_R = matrix(state[((2*num_patch)+1):(3*num_patch)], ncol = 1)
        P_S = matrix(state[((3*num_patch)+1):(4*num_patch)], ncol = 1)
        P_I = matrix(state[((4*num_patch)+1):(5*num_patch)], ncol = 1)
        S_S = matrix(state[((5*num_patch)+1):(6*num_patch)], ncol = 1)
        S_I = matrix(state[((6*num_patch)+1):(7*num_patch)], ncol = 1)
        
        ###Demographic parameters
        b_H <- matrix(rep(param["b_H"],num_patch), ncol = 1) #Human birth rate
        b_P <- matrix(rep(param["b_P"],num_patch), ncol = 1)  #P.vector birth rate
        b_S <-  matrix(rep(param["b_S"],num_patch), ncol = 1)  #S. vector birth rate
        mu_H <- matrix(rep(param["mu_H"],num_patch), ncol = 1)  #Human death rate
        mu_P <-  matrix(rep(param["mu_P"],num_patch), ncol = 1)  #P. vector death rate
        mu_S <- matrix(rep(param["mu_S"],num_patch), ncol = 1)  #S. vector death rate
        
        #Force of infection parameters
        a_P <- matrix(rep(param["a_P"],num_patch) , ncol = 1)  #biting rate of the p. vector
        a_S <-  matrix(rep(param["a_S"],num_patch) , ncol = 1)  #biting rate of the s.vector
        
        phi_P <-  matrix(rep(param["phi_P"],num_patch) , ncol = 1)  #transmission probability of p. vector
        phi_S <-  matrix(rep(param["phi_S"],num_patch) , ncol = 1)   #transmission probability of s. vector
        phi_H  <-  matrix(rep(param["phi_H"],num_patch) , ncol = 1)   #transmission probability of human
        
        # Recovery rate
        gamma <-  matrix(rep(param["gamma"],num_patch), ncol = 1)   #recovery rate of infected human
        
        #competition coefficient
        c_PS <-  matrix(rep(param["c_PS"],num_patch), ncol = 1)   #competitition effect of p.vector on s.vector
        c_SP <-  matrix(rep(param["c_SP"],num_patch), ncol = 1)  #competition effect of s.vector on p.vector
        
        ### FOI
        FOI_P <- matrix(a_P * phi_P, ncol = 1) #FOI for a primary vector
        FOI_S <- matrix(a_S * phi_S , ncol = 1) #FOI for a secondary vector
        FOI_H_P <- matrix(a_P * phi_H,ncol =1) #FOI for a human to a primary vector
        FOI_H_S <- matrix(a_S * phi_S,ncol = 1) #FOI for a human to a secondary vector
        
        ###Dispersal matrix (density-dependent function)
        a_max <- param["a_max"]
        k <- param["k"]#steepness parameter
        a_0 <- param["a_0"]  ###midpoint population size

        ###Dispersal matrix
        ###exponential decay parameter
        lambda <- param["lambda"]
        
        
        ### Population size
        N_H <- matrix(H_S + H_I + H_R)#human host population 
        N_P <- matrix(P_S + P_I) #p. vector population
        N_S <- matrix(S_S + S_I) #s. vector population
        N_V <- matrix(N_P + N_S) # primary and secondary vector population
       
         ###Human host 
        ###Susceptible 
        dH_S <- b_H * (N_H) - (FOI_P * H_S * (P_I/N_P)) - (FOI_S *H_S*(S_I/N_S))- (mu_H*H_S) 
        
        ###Infected
        dH_I <- (FOI_P * H_S * (P_I/N_P)) + (FOI_S *H_S*(S_I/N_S))- (gamma*H_I) -(mu_H*H_I)
        
        ###Recovered
        dH_R <- (gamma*H_I)-(mu_H*H_R)
        
        ###Account for the distance...
        
        probability_matrix <- apply(adj_matrix, c(1, 2),function(x) {
                if (x != 0) {
                        return(exp(-lambda * x))
                } else {
                        return(x)
                }
        })
        
        
        ###Account for the density-dependence
        dd_mat <- diag(c(a_max/(1.00 +  exp(-k *(N_V - a_0 )))),num_patch,num_patch)
        adjusted_prob <- sweep(probability_matrix, 1, rowSums(probability_matrix), FUN = "/")
        adjusted_prob[is.nan(adjusted_prob )== TRUE] = 0
        adjusted_prob_2 <- adjusted_prob * probability_matrix
        disp.contact2 <-adjusted_prob_2 %*% dd_mat



        ###P. vector
        ###Susceptible
        dP_S <-  (b_P* N_P) - (FOI_H_P *P_S* (H_I/N_H)) - (mu_P *P_S) - c_SP*(P_S)*(N_S) -
                (colSums(disp.contact2) * P_S) +  (disp.contact2 %*% P_S)
                

        ###Infected
        dP_I <- (FOI_H_P * P_S * (H_I/N_H)) - mu_P*P_I- c_SP*(P_I)*(N_S) #-
                (colSums(disp.contact2 ) * P_I) +  (disp.contact2 %*% P_I)
        
        
        ###S. vector 
        ###Susceptible
        dS_S <-  (b_S * N_S) - (FOI_H_S *S_S* (H_I/N_H)) - mu_S *S_S -  c_PS*(S_S)*(N_P) -
                (colSums(disp.contact2 ) * S_S) +  (disp.contact2 %*% S_S)
        
        
        ###Infected
        dS_I <- (FOI_H_S * S_S * (H_I/N_H)) - mu_S*S_I -  c_PS*(S_I)*(N_P) -
               (colSums(disp.contact2) * S_I) +  (disp.contact2 %*% S_I)
        
        
        return(list(c(dH_S,dH_I,dH_R,dP_S, dP_I, dS_S, dS_I)))
        }
        )
}

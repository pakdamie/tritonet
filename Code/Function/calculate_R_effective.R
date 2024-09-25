calculate_R_effective <- function() 
        
mu_H = 1/27375 #Human death rate
mu_P = 0.05 #P. vector death rate
mu_S = 0.05 #S. vector death rate

a_P = 7 #daily biting rate of the p. vector
a_S = 7 * (0.55) #daily biting rate of the s.vector

phi_P = 0.132 #transmission probability of p. vector
phi_S = 0.132 * 0.55 #transmission probability of s. vector
phi_H  = 0.116 #transmission probability of human

# Recovery rate of the acute phase
gamma = 1/56  #recovery rate of infected human

#competition coefficient
c_PS = 5e-4 #competitition effect of p.vector on s.vector
c_SP = 1e-6#competitition effect of s.vector on p.vector

a_max = 2 * 10^-2 
k = 0.02
a_0 = 300

H_I <- sample((1:10), 2, replace = TRUE)
H_S <- sample((1:10), 2, replace = TRUE)
H_R <- sample((1:10), 2, replace = TRUE)



P_S <- sample((1:10), 2, replace = TRUE)
P_I <- sample((1:10), 2, replace = TRUE)
S_S <- sample((1:10), 2, replace = TRUE)
S_I <- sample((1:10), 2, replace = TRUE)

N_S <- S_I + S_S
N_P <- P_I + P_S
N_H <- H_S + H_I + H_R

adjacency_matrix <- matrix(c(0,1,1,0), 
                           ncol = 2, 
                           byrow = TRUE)

Dispersal_Matrix <- a_max/(1 + k*(( (N_P + N_S) -a_0)))


HI_Input_from_P <- ((phi_P * a_P) * P_I[1]/N_P[1]) 
HI_Input_from_S <-  ((phi_S* a_S) * S_I[1]/N_S[1]) 

PI_Input_from_H <- ((phi_H * a_P) * H_I[1]/N_H[1]) + (adjacency_matrix %*%((Dispersal_Matrix) * matrix(P_I ,ncol = 1)))[1]
SI_Input_from_H <- ((phi_H * a_S) * H_I[1]/N_H[1]) + (adjacency_matrix %*%((Dispersal_Matrix) * matrix(S_I,ncol = 1)))[1]
          
                                                        
                                                                                               
HI_Output <- -(gamma + mu_H)
PI_Output <- -(mu_P + (c_SP*N_S[1])) - (adjacency_matrix %*%((Dispersal_Matrix) * matrix(P_I ,ncol = 1)))[1]
SI_Output <- -(mu_S + (c_PS *N_P[1]))   -(adjacency_matrix %*%((Dispersal_Matrix) * matrix(P_I ,ncol = 1)))[1]

F_matrix <- matrix(c(0,HI_Input_from_P ,HI_Input_from_S,
              PI_Input_from_H,0,0,
              SI_Input_from_H,0,0),ncol = 3, byrow = TRUE)
              
V_matrix <- matrix(c(HI_Output,0,0,
                     0,PI_Output,0,
                     0,0,SI_Output), ncol =3,byrow = TRUE)

max(eigen(F_matrix%*%solve(V_matrix))$values)

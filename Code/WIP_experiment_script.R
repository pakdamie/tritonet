###List of the 
low_connectance_network <- simulate_spatial_network(24601, 20, 0.05)
mid_connectance_network <- simulate_spatial_network(24601, 20, 0.10)

param_ALT<- c(b_H = (1/(300)), #Human mortality rate
           b_P = 0.15, # P. Vector birth rate
           b_S = 0.15,  # S. Vector birth rate
           mu_H = (1/(300)),  #Human death rate
           mu_V = 0,  #P. vector death rate
           f_P = 0.50, # Biting rate of the p. vector
           f_S = 0.25, # Biting rate of the s.vector
           theta_P = 0.40, # Transmission probability of p. vector
           theta_S = 0.20, # Transmission probability of s. vector
           theta_H  = 0.50, # Transmission probability of human
           gamma = 1/30,  #Recovery rate of infected human
           c_PS = 7e-5, #Competition effect of p.vector on s.vector
           c_SP = 1e-5,  #Competition effect of s.vector on p.vector
           c_PP = 1e-4, #Competition effect of p.vector on s.vector
           c_SS = 9e-5,
           d = 0.10
)



test_function_L_1 <- discrete_trito_model(1000, param_ALT,
                                           500, 0.25, 0.95 ,0.50,
                                           low_connectance_network, 1)

test_function_L_2 <-discrete_trito_model(1000,param_ALT,
                                        500, 0.25, 0.95,0.50,
                                       mid_connectance_network, 1)


plot_list_groups(test_function_L_1)[[1]] + plot_list_groups(test_function_L_2)[[1]]


calculate_R_effective_discrete_net(param_ALT,test_function_L_1,500)
calculate_R_effective_discrete_patch(param_ALT,test_function_L_1,500)
calculate_R_effective_discrete_net(param_ALT,test_function_L_2,500)
calculate_R_effective_discrete_patch(param_ALT,test_function_L_2,500)


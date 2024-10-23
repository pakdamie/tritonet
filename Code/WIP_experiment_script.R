###List of the 
low_connectance_network <- simulate_spatial_network(24601, 20, 0.05)
mid_connectance_network <- simulate_spatial_network(24601, 20, 0.10)

param <- c(b_H = (1/(100)),
           b_P =0.25, #P.vector birth rate
           b_S = 0.25,  #S. vector birth rate
           mu_H = (1/(100)),  #Human death rate
           mu_V = 0,  #P. vector death rate
           f_P = 0.50, #daily biting rate of the p. vector
           f_S = 0.25, #daily biting rate of the s.vector
           theta_P = 0.70, #transmission probability of p. vector
           theta_S = 0.50, #transmission probability of s. vector
           theta_H  = 0.90, #transmission probability of human
           # Recovery rate of the acute phase
           gamma = 1/50,  #recovery rate of infected human
           #competition coefficient
           c_PS = 3e-5,#competition effect of p.vector on s.vector
           c_SP = 1e-6,  #competition effect of s.vector on p.vector
           c_PP  = 5e-5,
           c_SS = 2e-5,
           d = 0.10
)



test_function_L_1 <- discrete_trito_model(1000, param,
                                           500, 0.45, 0.9 , 0.25,
                                           low_connectance_network, 1)


test_function_L_2 <-discrete_trito_model(1000, param,
                                       500, 0.4, 0.9,0.25,
                                       mid_connectance_network, 1)


plot_list_groups(test_function_L_1)[[1]] + plot_list_groups(test_function_L_2)[[1]]


###WORKING


param <- c(b_H = (1/(100)),
           b_P = (1), #P.vector birth rate
           b_S = (1),  #S. vector birth rate
           mu_H =  (1/(100)),  #Human death rate
           mu_V = 0.5,  #P. vector death rate
           f_P = 2.5, #daily biting rate of the p. vector
           f_S = 1.0, #daily biting rate of the s.vector
           theta_P = 0.25, #transmission probability of p. vector
           theta_S = 0.10, #transmission probability of s. vector
           theta_H  = 0.50, #transmission probability of human
           # Recovery rate of the acute phase
           gamma = 1/56,  #recovery rate of infected human
           #competition coefficient
           c_PS = 3e-5,#competition effect of p.vector on s.vector
           c_SP = 1e-6,  #competition effect of s.vector on p.vector
           c_PP  = 5e-5,
           c_SS = 5e-5,
           d = 0.5
)


test_function_L_1 <- discrete_trito_model(1000, param,
                                          500, 0.65, 0.9 ,0.50,
                                          low_connectance_network, 1)


plot(test_function_L_1[[1]][[4]][,3])

     
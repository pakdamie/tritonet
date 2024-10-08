

low_connectance_network <-  simulate_spatial_network(24601,30,0.1)
high_connectance_network <- simulate_spatial_network(24601,30,0.50)



param <- c(b_H = (1/(300)),
           b_P = (1.0), #P.vector birth rate
           b_S = (1.0),  #S. vector birth rate
           mu_H =  (1/(300)),  #Human death rate
           mu_V = 0.25,  #P. vector death rate
           f_P = 1.5, #daily biting rate of the p. vector
           f_S = 0.5, #daily biting rate of the s.vector
           theta_P = 0.25, #transmission probability of p. vector
           theta_S = 0.20, #transmission probability of s. vector
           theta_H  = 0.25, #transmission probability of human
           # Recovery rate of the acute phase
           gamma = 1/80,  #recovery rate of infected human
           #competition coefficient
           c_PS = 8e-5,#competition effect of p.vector on s.vector
           c_SP = 5e-6,  #competition effect of s.vector on p.vector
           c_PP  = 1e-4,
           c_SS = 9e-5,
           d = 0.15
)




test_function_L_1 <-discrete_trito_model(1000, param,
                      500, 0.65, 0.9 ,1,
                      low_connectance_network, .1)
test_function_L_2 <-discrete_trito_model(1000, param,
                                       500, 0.25, 0.9 ,1,
                                       low_connectance_network, .1)



plot_list_groups(list = test_function_L_1[[1]][[1]]) + ggtitle("HS - L")
plot_list_groups(list = test_function_L_1[[1]][[2]]) + ggtitle("HI - L")
plot_list_groups(list = test_function_L_1[[1]][[3]])
plot_list_groups(list = test_function_L_1[[1]][[4]]) + ggtitle("PS - L")
plot_list_groups(list = test_function_L_1[[1]][[5]]) + ggtitle("PI - L")
plot_list_groups(list = test_function_L_1[[1]][[6]]) +  ggtitle("SS - L")
plot_list_groups(list = test_function_L_1[[1]][[7]]) +  ggtitle("SI - L")


plot_list_groups(list = test_function_L[[1]][[4]] +  test_function_L[[1]][[5]]) 
plot_list_groups(list =  test_function_L[[1]][[6]] + test_function_L [[1]][[7]]) #SS


test_function_H <-discrete_trito_model(1000, param,
                                       500, 0.25, 0.9 ,0.25,
                                       high_connectance_network, .1)

plot_list_groups(list = test_function_H[[1]][[1]])
plot_list_groups(list = test_function_H[[1]][[2]])
plot_list_groups(list = test_function_H[[1]][[3]])
plot_list_groups(list = test_function_H[[1]][[4]]) #PS
plot_list_groups(list = test_function_H[[1]][[5]]) #PI
plot_list_groups(list =  test_function_H[[1]][[6]]) #SS
plot_list_groups(list = test_function_H [[1]][[7]]) #SI


plot_R0_groups(test_function_L_1[[2]]) + ylim(0.5,2) + 
plot_R0_groups(test_function_L_2[[2]]) + ylim(0.5,2)

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




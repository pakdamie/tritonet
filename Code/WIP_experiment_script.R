param <- c(b_H = (1/(50)),
           b_P = (1.2 + 1e-7), #P.vector birth rate
           b_S = 1.2,  #S. vector birth rate
           mu_H =  (1/(50)),  #Human death rate
           mu_V = 1.2,  #P. vector death rate
           f_P = 2.5, #daily biting rate of the p. vector
           f_S = 1.0, #daily biting rate of the s.vector
           theta_P = 0.05, #transmission probability of p. vector
           theta_S = 0.01, #transmission probability of s. vector
           theta_H  = 0.50, #transmission probability of human
           # Recovery rate of the acute phase
           gamma = 1/56,  #recovery rate of infected human
           #competition coefficient
           c_PS = 1e-4,#competition effect of p.vector on s.vector
           c_SP = 1e-7,  #competition effect of s.vector on p.vector
           d = 0.25,
           k_p = 1500,
           k_s = 1500
)


low_connectance_network <-  simulate_spatial_network(24601,30,0.01)
high_connectance_network <- simulate_spatial_network(24601,30,0.50)



test_function_L <-discrete_trito_model(200, param,
                        100, 0.4, 0.7 ,0.25,
                      low_connectance_network )

test_function_H <-discrete_trito_model(200, param,
                                       100, 0.4, 0.7 ,0.25,
                                      high_connectance_network )


plot_list_groups(list = test_function_L[[1]][[1]]) + ggtitle("HS - L")
plot_list_groups(list = test_function_L[[1]][[2]]) + ggtitle("HI - L")
plot_list_groups(list = test_function_L[[1]][[3]])
plot_list_groups(list = test_function_L[[1]][[4]]) #PS
plot_list_groups(list = test_function_L[[1]][[5]]) #PI
plot_list_groups(list =  test_function_L[[1]][[6]]) #SS
plot_list_groups(list = test_function_L [[1]][[7]]) #SI

plot_list_groups(list = test_function_H[[1]][[1]])
plot_list_groups(list = test_function_H[[1]][[2]])
plot_list_groups(list = test_function_H[[1]][[3]])
plot_list_groups(list = test_function_H[[1]][[4]]) #PS
plot_list_groups(list = test_function_H[[1]][[5]]) #PI
plot_list_groups(list =  test_function_H[[1]][[6]]) #SS
plot_list_groups(list = test_function_H [[1]][[7]]) #SI


plot_R0_groups(test_function_L[[2]])
plot_R0_groups(test_function_H[[2]])




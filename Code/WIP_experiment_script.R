low_connectance_network <- simulate_spatial_network(24601, 20, 0.1)


param <- c(b_H = (1/(100)),
           b_P = (1.10), #P.vector birth rate
           b_S = (1.10),  #S. vector birth rate
           mu_H =  (1/(100)),  #Human death rate
           mu_V = 0,  #P. vector death rate
           f_P = 1.5, #daily biting rate of the p. vector
           f_S = 1.0, #daily biting rate of the s.vector
           theta_P = 0.15, #transmission probability of p. vector
           theta_S = 0.10, #transmission probability of s. vector
           theta_H  = 0.20, #transmission probability of human
           # Recovery rate of the acute phase
           gamma = 1/100,  #recovery rate of infected human
           #competition coefficient
           c_PS = 8e-5,#competition effect of p.vector on s.vector
           c_SP = 1e-5,  #competition effect of s.vector on p.vector
           c_PP  = 1e-4,
           c_SS = 9e-5,
           d = 0.50
)


test_function_L_1<-discrete_trito_model(1000, param,
                                           500, 0.65, 0.9 ,0.50,
                                           low_connectance_network, 0.1)

test_function_L_2 <-discrete_trito_model(1000, param,
                                       500, 0.25, 0.9 ,0.50,
                                       low_connectance_network, 0.1)
test_function_L_3 <-discrete_trito_model(1000, param,
                                         5000, 0.25, 0.9 ,0.50,
                                         low_connectance_network, 0.1)


plot_list_groups(list = test_function_L_2[[1]][[1]]) + ggtitle("HS - L") +
plot_list_groups(list = test_function_L_2[[1]][[2]]) + ggtitle("HI - L") + 
plot_list_groups(list = test_function_L_2[[1]][[3]]) + ggtitle("HR - L") + 
plot_list_groups(list = test_function_L_2[[1]][[4]]) + ggtitle("PS - L") + 
plot_list_groups(list = test_function_L_2[[1]][[5]]) + ggtitle("PI - L") + 
plot_list_groups(list = test_function_L_2[[1]][[6]]) +  ggtitle("SS - L") + 
plot_list_groups(list = test_function_L_2[[1]][[7]]) +  ggtitle("SI - L") 



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


plot_R0_groups(test_function_L_1[[2]]) + ylim(0.5,2) + ggtitle("65% of primary survive") + 
plot_R0_groups(test_function_L_2[[2]]) + ylim(0.5,2) + ggtitle("25% of primary survive") + 
plot_R0_groups(test_function_L_3[[2]]) + ylim(0.5,2) + ggtitle("No disturbance")


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



c_PS = 8e-5,#competition effect of p.vector on s.vector
c_SP = 5e-6,  #competition effect of s.vector on p.vector
c_PP  = 1e-4,
c_SS = 9e-5,



connectance_010_network <- simulate_spatial_network(24601, 20, 0.1)
connectance_015_network <- simulate_spatial_network(24601, 20, 0.15)
connectance_020_network <- simulate_spatial_network(24601, 20, 0.2)
connectance_025_network <- simulate_spatial_network(24601, 20, 0.25)
connectance_030_network <- simulate_spatial_network(24601, 20, 0.3)

test_function_010<-discrete_trito_model(1000, param,
                                        500, 0.25, 0.9, 0.50,
                                        connectance_010_network , 0.5)

test_function_020<-discrete_trito_model(1000, param,
                                        500, 0.25, 0.8 ,0.50,
                                        connectance_020_network ,0.5)

test_function_030<-discrete_trito_model(1000, param,
                                        500, 0.25, 0.8 ,0.50,
                                        connectance_030_network , 0.5)
test_function_040<-discrete_trito_model(1000, param,
                                        500, 0.25, 0.8 ,0.50,
                                        connectance_040_network ,0.5)
test_function_050 <-discrete_trito_model(1000, param,
                                         500, 0.25, 0.8 ,0.50,
                                         connectance_050_network , 0.5)


R0_max_010 <- cbind(calculate_max_CV_RE(test_function_010[[2]]), id = 0.10) 
R0_max_020 <- cbind(calculate_max_CV_RE(test_function_020[[2]]), id = 0.20) 
R0_max_030 <- cbind(calculate_max_CV_RE(test_function_030[[2]]), id = 0.30) 
R0_max_040 <- cbind(calculate_max_CV_RE(test_function_040[[2]]), id = 0.40) 
R0_max_050 <- cbind(calculate_max_CV_RE(test_function_050[[2]]), id = 0.50) 

R0_max_all <- rbind(
        R0_max_010,
        R0_max_020,
        R0_max_030,
        R0_max_040,
        R0_max_050)

Max_R_Effective <- ggplot(R0_max_all, aes( x= id, y = max)) + geom_line(size = 1.2) + xlab("Connectance") + 
        ylab("Maximum R-eff") + theme_classic() + 
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 15)) + 
        ggtitle("Maximum R effective") 

CV_R_Effective <- ggplot(R0_max_all, aes( x= id, y = cv)) + geom_line(size = 1.2) + xlab("Connectance") + 
        ylab("CV R-eff") + theme_classic() + 
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 15)) + 
        ggtitle("CV R effective") 

Max_R_Effective / CV_R_Effective


param_alt <- c(
  b_H = 1 / (365 * 70), ## Human mortality rate
  b_P = 1/10, # P. Vector birth rate
  b_M = 1/10, # S. Vector birth rate
  mu_H = 1 / (365 * 70), ## Human death rate
  f_P = 0.25, # Biting rate of the p. vector
  f_M = 0.25 * 0.75, # Biting rate of the s.vector
  theta_P = 0.50, # Transmission probability of p. vector
  theta_M = 0.50 * 0.75, # Transmission probability of s. vector
  theta_H = 0.50, # Transmission probability of human
  gamma = 1 / 90, # Recovery rate of infected human
  c_PM = 3e-4, ## Competition effect of p.vector on s.vector
  c_MP = 3e-6, ## Competition effect of s.vector on p.vector
  c_PP = 4.5e-4, ## Competition effect of p.vector on s.vector
  c_MM = 2.5e-4, ## Competition effect of s.vector on s.vector
  mu_V = 1/20,
  ntime = (365 * 50),
  disturbance_time = (365 * 25),
  delta_T = 1,
  prop = 1,
  mortality_P = 0.50, # This will change
  mortality_M = 0.9
)


Initial_List <- create_initial_states(param_alt)

model_output_list <-
  discrete_trito_model_rcpp_ONEPATCH_alt(
    HS = Initial_List[[1]],
    HI = Initial_List[[2]],
    HR = Initial_List[[3]],
    PS = Initial_List[[4]],
    PI = Initial_List[[5]],
    MS = Initial_List[[6]],
    MI = Initial_List[[7]],
    param = param_alt
  ) 

RE_output <- Calculate_Human_REff(model_output_list,param = param_alt)


plot_list_groups(model_output_list)

r#Plot total vector abundance over time with the RE As color

ggplot(RE_output  , aes(x = time, y= NP )) + geom_path(color = 'red') + 
  geom_path(aes( x= time, y = NM), color = 'blue')

a<- (ggplot(subset(RE_output,
              RE_output$time > 9100 &
              RE_output$time < 9300),aes(x = time, y= Primary_HP))  +
  geom_path(color = 'red')) + 
  geom_path(data=subset(RE_output,
                   RE_output$time > 9100 &
                     RE_output$time < 9300),
            aes(x = time, y= Secondary_HM)) + 
  
 b<- (ggplot(subset(RE_output,
                RE_output$time > 9100 &
                  RE_output$time < 9300),aes(x = time, y= Primary_PH)) + 
  geom_path(color = 'red') ) +
  geom_path(data=subset(RE_output,
                        RE_output$time > 9100 &
                          RE_output$time < 9300),
            aes(x = time, y= Secondary_MH))
a/b

c <-  (ggplot(subset(RE_output,
                     RE_output$time > 9100 &
                       RE_output$time < 9300),aes(x = time, y=  PtoH/RE)) + 
         geom_path(color = 'red') ) +
  geom_path(data=subset(RE_output,
                        RE_output$time > 9100 &
                          RE_output$time < 9300),
            aes(x = time, y=  MtoH/RE))
(ggplot(subset(RE_output,
               RE_output$time > 9100 &
                 RE_output$time < 9300),aes(x = time, y=  HS)) + 
    geom_path(color = 'red') ) + 
  geom_path()

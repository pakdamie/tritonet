# This is an analysis script for the one patch dynamics.

param_same <-
  c(
    b_H = 1/ (6570), ## Human mortality rate
    b_P = 0.1, # P. Vector birth rate
    b_M = 0.1, # S. Vector birth rate
    mu_H = 1/ (6570), ## Human death rate
    f_P = 0.05, # Biting rate of the p. vector
    f_M = 0.05, # Biting rate of the s.vector
    theta_P = 0.90, # Transmission probability of p. vector
    theta_M = 0.90, # Transmission probability of s. vector
    theta_H = 0.50, # Transmission probability of human
    gamma = 1 / 100, # Recovery rate of infected human
    c_PM = 2e-5, ## Competition effect of p.vector on s.vector
    c_MP = 1e-5, ## Competition effect of s.vector on p.vector
    c_PP = 3e-5, ## Competition effect of p.vector on s.vector
    c_MM = 3e-5, ## Competition effect of s.vector on s.vector
    ntime = 365 * 50,
    disturbance_time = 365 * 25,
    delta_T = 1,
    d = 0, # Dispersal is automatically turned off, but just to make sure
    prop = 1,
    mortality_P = 0.01, # This will change
    mortality_M = 1
  )


param_standard <-
  c(
    b_H = 1/ (1000), ## Human mortality rate
    b_P = 0.1, # P. Vector birth rate
    b_M = 0.1, # S. Vector birth rate
    mu_H = 1/ (1000), ## Human death rate
    f_P = 0.05, # Biting rate of the p. vector
    f_M = 0.030, # Biting rate of the s.vector
    theta_P = 0.90, # Transmission probability of p. vector
    theta_M = 0.75, # Transmission probability of s. vector
    theta_H = 0.50, # Transmission probability of human
    gamma = 1 / 100, # Recovery rate of infected human
    c_PM = 2e-5, ## Competition effect of p.vector on s.vector
    c_MP = 1e-5, ## Competition effect of s.vector on p.vector
    c_PP = 3e-5, ## Competition effect of p.vector on s.vector
    c_MM = 3e-5, ## Competition effect of s.vector on s.vector
    ntime = 365 * 50,
    disturbance_time = 365 * 25,
    delta_T = 1,
    d = 0, # Dispersal is automatically turned off, but just to make sure
    prop = 1,
    mortality_P = 0.01, # This will change
    mortality_M = 1
  ) # This will not change

### What if M is a better vector BUT it is being repressed?
param_M_better <-
  c(
    b_H = 1 / (1000), ## Human mortality rate
    b_P = 0.1, # P. Vector birth rate
    b_M = 0.1, # S. Vector birth rate
    mu_H = 1 / (1000), ## Human death rate
    f_P = 0.03, # Biting rate of the p. vector
    f_M = 0.05, # Biting rate of the s. vector
    theta_P = 0.75, # Transmission probability of p. vector
    theta_M = 0.90, # Transmission probability of s. vector
    theta_H = 0.50, # Transmission probability of human
    gamma = 1 / 100, # Recovery rate of infected human
    c_PM = 2e-5, ## Competition effect of p.vector on s.vector
    c_MP = 1e-5, ## Competition effect of s.vector on p.vector
    c_PP = 3e-5, ## Competition effect of p.vector on s.vector
    c_MM = 3e-5, ## Competition effect of s.vector on s.vector
    ntime = 365 * 50,
    disturbance_time = 365 * 25,
    delta_T = 1,
    d = 0, # Dispersal is automatically turned off, but just to make sure
    prop = 1,
    mortality_P = 0.01, # This will change
    mortality_M = 1
  ) # This will not change

# The survivorship of the primary vector
mortality_P <- seq(0.1, 1, 0.10)

### Creating a list of parameters that I can loop through
parameter_mortality_P_list <- NULL
parameter_mortality_P_list_switch <- NULL
parameter_mortality_P_list_same <- NULL

for (i in seq(1, length(mortality_P))) {
  copied_param <- param_standard
  copied_param["mortality_P"] <- mortality_P[i]
  ## Put in list to be looped through
  parameter_mortality_P_list[[i]] <- copied_param

  copied_param_switch <- param_M_better
  copied_param_switch["mortality_P"] <- mortality_P[i]
  ## Put in list to be looped through
  parameter_mortality_P_list_switch[[i]] <- copied_param_switch
  
  copied_param_same <- param_same
  copied_param_same["mortality_P"] <- mortality_P[i]
  ## Put in list to be looped through
  parameter_mortality_P_list_same[[i]] <- copied_param_same
  
  
}

### Setting initial conditions
Initial_List <- create_initial_states(param_standard, patch_num = 1)
Initial_List_NoSecondary <- Initial_List
Initial_List_NoSecondary[[6]][1] <- 0
Initial_List_NoSecondary[[7]][1] <- 0


### Simulate the model with varying primary-removal
model_output_mort_P <- NULL
model_output_mort_P_switch <- NULL
model_output_mort_P_noM <- NULL
model_output_mort_P_same <- NULL
for (p in seq(1, length(parameter_mortality_P_list))) {
  model_output_mort_P[[p]] <-
    discrete_trito_model_rcpp_ONEPATCH(
      HS = Initial_List[[1]],
      HI = Initial_List[[2]],
      HR = Initial_List[[3]],
      PS = Initial_List[[4]],
      PI = Initial_List[[5]],
      MS = Initial_List[[6]],
      MI = Initial_List[[7]],
      param = parameter_mortality_P_list[[p]]
    )

  model_output_mort_P_switch[[p]] <-
    discrete_trito_model_rcpp_ONEPATCH(
      HS = Initial_List[[1]],
      HI = Initial_List[[2]],
      HR = Initial_List[[3]],
      PS = Initial_List[[4]],
      PI = Initial_List[[5]],
      MS = Initial_List[[6]],
      MI = Initial_List[[7]],
      param = parameter_mortality_P_list_switch[[p]]
    )

  model_output_mort_P_noM[[p]] <-
    discrete_trito_model_rcpp_ONEPATCH(
      HS = Initial_List_NoSecondary[[1]],
      HI = Initial_List_NoSecondary[[2]],
      HR = Initial_List_NoSecondary[[3]],
      PS = Initial_List_NoSecondary[[4]],
      PI = Initial_List_NoSecondary[[5]],
      MS = Initial_List_NoSecondary[[6]],
      MI = Initial_List_NoSecondary[[7]],
      param = parameter_mortality_P_list[[p]]
    )
  
  model_output_mort_P_same[[p]] <-
    discrete_trito_model_rcpp_ONEPATCH(
      HS = Initial_List[[1]],
      HI = Initial_List[[2]],
      HR = Initial_List[[3]],
      PS = Initial_List[[4]],
      PI = Initial_List[[5]],
      MS = Initial_List[[6]],
      MI = Initial_List[[7]],
      param = parameter_mortality_P_list_same[[p]]
    )
}


plot_list_groups(model_output_mort_P[[5]][1:7])
plot_list_groups(model_output_mort_P_noM[[1]][1:7])












## Calculate the RE over time

RE_onepatch_List <- NULL
RE_onepatch_List_Switch <- NULL
RE_onepatch_List_noM <- NULL
RE_onepatch_List_Same <- NULL


for (p in seq(1, length(parameter_mortality_P_list))) {
 
  RE_onepatch_List[[p]] <- cbind.data.frame(
    calculate_RE_onepatch(
      model_output_mort_P[[p]], # Model outputs
      parameter_mortality_P_list[[p]],
      1
    ), # Parameter choices
    primary_removal = mortality_P[[p]]
  ) # Factor
  
  RE_onepatch_List_Switch[[p]] <- cbind.data.frame(
    calculate_RE_onepatch(
      model_output_mort_P_switch[[p]], # Model outputs
      parameter_mortality_P_list_switch[[p]],
    1), # Parameter choices
    primary_removal = mortality_P[[p]]
  )
  
  RE_onepatch_List_noM[[p]] <- cbind.data.frame(
    calculate_RE_onepatch(
      model_output_mort_P_noM[[p]], # Model outputs
      parameter_mortality_P_list[[p]],
      0
    ), # Parameter choices
     primary_removal = mortality_P[[p]]
  )

  
  RE_onepatch_List_Same[[p]] <- cbind.data.frame(
    calculate_RE_onepatch(
      model_output_mort_P_same[[p]], # Model outputs
      parameter_mortality_P_list_same[[p]],
    1), # Parameter choices
    primary_removal = mortality_P[[p]]
  )
  
}

RE_onepatch_DF <- do.call(rbind, RE_onepatch_List)
RE_onepatch_DF_Switch <- do.call(rbind, RE_onepatch_List_Switch)
RE_onepatch_DF_noM <- do.call(rbind, RE_onepatch_List_noM)
RE_onepatch_DF_Same <- do.call(rbind, RE_onepatch_List_Same)

colnames(RE_onepatch_DF) <- c("RE", "N_P", "N_M", "H_S","P_S","M_S", 
                              "H_I","P_I","M_I",
                              "time","primary_removal")
colnames(RE_onepatch_DF_Switch) <- c("RE", "N_P", "N_M", "H_S","P_S","M_S", 
                                     "H_I","P_I","M_I",
                                     "time","primary_removal")
colnames(RE_onepatch_DF_noM) <-c("RE", "N_P", "N_M", "H_S","P_S","M_S", 
                                 "H_I","P_I","M_I",
                                 "time","primary_removal")
colnames(RE_onepatch_DF_Same) <- c("RE", "N_P", "N_M", "H_S","P_S","M_S", 
                                   "H_I","P_I","M_I",
                                   "time","primary_removal")


### What is the initial abundance at equilibrium?
# RE_onepatch_DF[RE_onepatch_DF$time ==  (365 * 25)- 1,]$N_P  #2857.143
# RE_onepatch_DF[RE_onepatch_DF$time ==  (365 * 25)- 1,]$N_M  #1428.571
# Full abundance: 2857.143 + 1428.571 = 4285.714
# Figure showing how the RE changes over time with
# changing total secondary vectors primary vectors

standard_PvsM_GG <-
  ggplot() +
  geom_path(data = 
    subset(
    RE_onepatch_DF,
    RE_onepatch_DF$time > 9123),
    aes(
      x = N_P, y = N_M, color = RE,
      group = primary_removal),size = 1.1) +
  geom_point(
    data =
      subset(
        RE_onepatch_DF,
        RE_onepatch_DF$time == 5000
      ),
    aes(x = N_P, y = N_M),
    size = 2.5,
    color = "black"
  ) +
  scale_colour_gradient2(
    low = "#007dfa",
    mid = "#f0f5f7",
    high = "#ff0d55",
    midpoint = 1,
    name = expression(R[E]),
    limit = c(0.8,1.25)
    
  ) +
  scale_x_continuous(breaks = seq(0,4000,500))+ 
  theme_classic() +
  xlab("Primary vector abundance") +
  ylab("Secondary vector abundance") +
  ggtitle("P.vector is the \nmore effective transmitter") + 
  theme(legend.position = "right") + 
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 13, color = 'black'))


ggsave(here("Figures_Process","Standard_Primary_Secondary_Traj_RE.pdf"), units = 'in',
       width = 6, height =4.5)



switch_PvsM_GG <-
  ggplot(
  subset(
    RE_onepatch_DF_Switch,
    RE_onepatch_DF_Switch$time > 9123
  ),
  aes(
    x = N_P, y = N_M, color = RE,
    group = primary_removal
  )
) +
  geom_path(size = 1.1) +
  geom_point(
    data =
      subset(
        RE_onepatch_DF_Switch,
        RE_onepatch_DF_Switch$time == 5000
      ),
    aes(x = N_P, y = N_M),
    size = 2.5,
    color = "black"
  ) +
  scale_colour_gradient2(
    low = "#007dfa",
    mid = "#f0f5f7",
    high = "#ff0d55",
    midpoint = 1,
    name = expression(R[E]),
    limit = c(0.8,1.25)
    
  ) +
  scale_x_continuous(breaks = seq(0,4000,500))+ 
  theme_classic() +
  xlab("Primary vector abundance") +
  ylab("Secondary vector abundance") +
  ggtitle("S.vector is the \nmore effective transmitter") + 
  theme(legend.position = "right") + 
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 13, color = 'black'))

same_PvsM_GG <- 
  ggplot(
  subset(
    RE_onepatch_DF_Same,
    RE_onepatch_DF_Same$time > 9123
  ),
  aes(
    x = N_P, y = N_M, color = RE,
    group = primary_removal
  )
) +
  geom_path(size = 1.1) +
  geom_point(
    data =
      subset(
        RE_onepatch_DF_Same,
        RE_onepatch_DF_Same$time == 5000
      ),
    aes(x = N_P, y = N_M),
    size = 2.5,
    color = "black"
  ) +
  scale_colour_gradient2(
    low = "#007dfa",
    mid = "#f0f5f7",
    high = "#ff0d55",
    midpoint = 1,
    name = expression(R[E]),
    limit = c(0.8,1.25)
  ) +
  scale_x_continuous(breaks = seq(0,4000,500))+ 
  theme_classic() +
  xlab("Primary vector abundance") +
  ylab("Secondary vector abundance") +
  ggtitle("Both are effective transmitter") + 
  theme(legend.position = "right") + 
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 13, color = 'black'))

standard_PvsM_GG + switch_PvsM_GG + same_PvsM_GG + plot_layout(guides = 'collect')


ggsave(here("Figures_Process","Primary_Secondary_Traj_RE.pdf"), units = 'in',
            width = 12, height =3)

###What about the primary vector only?


### The equilibirum is at 3333 (when you don't have any secondary
###vector)
ggplot( subset(
  RE_onepatch_DF_noM,
  RE_onepatch_DF_noM$time > 9123 &
  RE_onepatch_DF_noM$time  < 9200),
  aes( x= time - 9124, y= RE, color = N_P, group  = primary_removal)) + 
  geom_path(size = 1) + 
  scale_color_viridis(name = "Primary vector\nabundance") + 
  xlab("Time since disturbance") + 
  ylab(expression(R[E])) +
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15)) + 


ggplot(subset(
  RE_onepatch_DF_noM,
  RE_onepatch_DF_noM$time > 9123 &
    RE_onepatch_DF_noM$time  < 9200),
  aes( x= time - 9124, y= N_P, color = RE, group  = primary_removal)) + 
  geom_path(size = 1) + 
  scale_colour_gradient2(
    low = "#007dfa",
    mid = "#f0f5f7",
    high = "#ff0d55",
    midpoint = 1,
    name = expression(R[E]),
    limit = c(1,1.05)
  )  +
  xlab("Time since disturbance") + 
  ylab("Primary vector abundance") +
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15))

ggsave(here("Figures_Process","Primary_only_plots.pdf"),
       units = 'in', height = 4, width = 10)

### Is it just the rise of susceptibles as Caitlin suggests?


ggplot(subset(
  RE_onepatch_DF_noM,
  RE_onepatch_DF_noM$time > 9123 &
    RE_onepatch_DF_noM$time  < 9200),
  aes( x= time - 9124, y= P_S, color = RE, group  = primary_removal)) + 
  geom_path(size = 1) + 
  scale_colour_gradient2(
    low = "#007dfa",
    mid = "#f0f5f7",
    high = "#ff0d55",
    midpoint = 1,
    name = expression(R[E]),
    limit = c(1,1.05)
  )  +
  xlab("Time since disturbance") + 
  ylab("Primary vector abundance") +
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15))

###Number of new infectives for the humans
RE_onepatch_DF_noM$new_human_cases = 
  param_standard["f_P"] * param_standard["theta_P"] *
  RE_onepatch_DF_noM$P_I * 
  (RE_onepatch_DF_noM$H_S/1000)

RE_onepatch_DF_noM$new_primary_cases = 
  param_standard["f_P"] * param_standard["theta_H"] *
  RE_onepatch_DF_noM$H_I * 
  (RE_onepatch_DF_noM$P_S/1000)

ggplot(subset(
  RE_onepatch_DF_noM,
  RE_onepatch_DF_noM$time > 9123 &
    RE_onepatch_DF_noM$time  < 9300),
  aes( x= time - 9124, y= new_human_cases, color = RE, group  = primary_removal)) + 
  geom_path(size = 1) +geom_path(size = 1) + 
  scale_colour_gradient2(
    low = "#007dfa",
    mid = "#f0f5f7",
    high = "#ff0d55",
    midpoint = 1,
    name = expression(R[E]),
    limit = c(1,1.05)
  )  +
  xlab("Time since disturbance") + 
  ylab("New human infection cases") +
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15))


ggplot(subset(
  RE_onepatch_DF_noM,
  RE_onepatch_DF_noM$time > 9123 &
    RE_onepatch_DF_noM$time  < 9300),
  aes( x= P_S, y= new_primary_cases, color = N_P, group  = primary_removal)) + 
  geom_path(size = 1) +
  geom_path(size = 1)  +
  xlab("New human infection cases") + 
  ylab("RE") +
  theme_classic() + 
  scale_color_viridis() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15)) + 
  facet_wrap(~primary_removal)


# RE_onepatch_DF_noM[RE_onepathsearch_db_keywords()# RE_onepatch_DF_noM[RE_onepatch_DF_noM$time ==  (365 * 25)- 1,]$N_P  #3333.333

###Are there any points where the abundance is greater than 4285.714
RE_onepatch_DF_subset <- subset(
  RE_onepatch_DF,
  RE_onepatch_DF$time > 9124 & 
    RE_onepatch_DF$N_P + RE_onepatch_DF$N_M > 4285.714)

RE_onepatch_DF_Switch_subset <- subset(
  RE_onepatch_DF_Switch,
  RE_onepatch_DF_Switch$time > 9124 & 
    RE_onepatch_DF_Switch$N_P + RE_onepatch_DF_Switch$N_M > 4285.714
)

RE_onepatch_DF_Same_subset <- subset(
  RE_onepatch_DF_Same,
  RE_onepatch_DF_Same$time > 9124 & 
  RE_onepatch_DF_Same$N_P + RE_onepatch_DF_Same$N_M > 4285.714
)
###




ggplot(subset(RE_onepatch_DF,
  RE_onepatch_DF$time > 9100  & RE_onepatch_DF$time < 9500 & 
    RE_onepatch_DF$primary_removal == 0.1 ),
       aes(x = time, y = RE, group = primary_removal, color = H_S)) + 
  geom_path(size = 1.1) + 
  scale_color_viridis(name = "Total vector abundance" ,option = 'viridis') + 
  theme_dark() + 
  xlab("Time") 

RE_onepatch_DF_01<- subset(RE_onepatch_DF,
       RE_onepatch_DF$time > 9100  & RE_onepatch_DF$time < 15000 & 
         RE_onepatch_DF$primary_removal == 0.1 )


RE_onepatch_DF_50<- subset(RE_onepatch_DF,
                           RE_onepatch_DF$time > 9100  & RE_onepatch_DF$time < 15000 & 
                             RE_onepatch_DF$primary_removal == 0.50 )


plot(RE_onepatch_DF_01$time,RE_onepatch_DF_01$H_S,
     type = 'l', xlab = "Time", ylab = 
       "Susceptible humans")
abline(v= 9125, col = 'red')


plot(RE_onepatch_DF_50$time,RE_onepatch_DF_50$H_S,
     type = 'l', xlab = "Time", ylab = 
       "Susceptible humans")
abline(v= 9125, col = 'red')

plot(RE_onepatch_DF_01$H_S,
     param_same['f_P'] * param_same['theta_P'] * RE_onepatch_DF_01$H_S/1000 *
       RE_onepatch_DF_01$P_I,type = 'l', xlab = "Susceptible humans",
     ylab = "New humans infections (P only)" )

plot(RE_onepatch_DF_50$H_S,
     param_same['f_P'] * param_same['theta_P'] * RE_onepatch_DF_50$H_S/1000 *
       RE_onepatch_DF_50$P_I,type = 'l', xlab = "Susceptible humans",
     ylab = "New humans infections (P only)" )

plot(RE_onepatch_DF_01$H_S,
     param_same['f_M'] * param_same['theta_M'] * RE_onepatch_DF_01$H_S/1000 *
       RE_onepatch_DF_01$M_I,type = 'l', xlab = "Susceptible humans",
     ylab = "New humans infections (M only)" )
plot(RE_onepatch_DF_50$H_S,
     param_same['f_M'] * param_same['theta_M'] * RE_onepatch_DF_50$H_S/1000 *
       RE_onepatch_DF_50$M_I,type = 'l', xlab = "Susceptible humans",
     ylab = "New humans infections (M only)" )

plot(RE_onepatch_DF_01$N_P + RE_onepatch_DF_01$N_M,
     (param_same['f_P'] * param_same['theta_P'] * RE_onepatch_DF_01$H_S/1000 *
      RE_onepatch_DF_01$P_I) + 
       (param_same['f_M'] * param_same['theta_M'] *RE_onepatch_DF_01$H_S/1000 *
       RE_onepatch_DF_01$M_I),type = 'l', xlab = "Total vector population",
     ylab = "New human infections (P AND M)")
abline(v= 9125, col = 'red')

###What if we switch it?


ggplot(subset(RE_onepatch_DF_Switch,
              RE_onepatch_DF_Switch$time > 9100  & 
                RE_onepatch_DF_Switch$time < 9500 & 
                RE_onepatch_DF_Switch$primary_removal == 0.1 ),
       aes(x = time, y = RE, group = primary_removal, color = N_P + N_M)) + 
  geom_path(size = 1.1) + 
  scale_color_viridis(name = "Total vector abundance" ,option = 'viridis') + 
  theme_dark() + 
  xlab("Time") 


RE_onepatch_Switch_DF_01<- subset(RE_onepatch_DF_Switch,
                                  RE_onepatch_DF_Switch$time > 9100  & 
                                    RE_onepatch_DF_Switch$time < 15000 & 
                                    RE_onepatch_DF_Switch$primary_removal == 0.1 )

plot(RE_onepatch_Switch_DF_01$H_S,
     param_same['f_M'] * param_same['theta_M'] * RE_onepatch_Switch_DF_01$H_S/1000 *
       RE_onepatch_Switch_DF_01$M_I,type = 'l', xlab = "Susceptible humans",
     ylab = "New humans infections (M only)" )

plot(RE_onepatch_Switch_DF_01$H_S,
     param_same['f_P'] * param_same['theta_P'] * RE_onepatch_Switch_DF_01$H_S/1000 *
       RE_onepatch_Switch_DF_01$P_I,type = 'l', xlab = "Susceptible humans",
     ylab = "New humans infections (P only)" )


ggplot(subset(RE_onepatch_DF_noM,
              RE_onepatch_DF_noM$time > 9100  & 
                RE_onepatch_DF_noM$time < 9500 & 
                RE_onepatch_DF_noM$primary_removal == 0.1 ),
       aes(x = time, y = RE, group = primary_removal, color = N_P + 0)) + 
  geom_path(size = 1.1) + 
  scale_color_viridis(name = "Total vector abundance" ,option = 'viridis') + 
  theme_dark() + 
  xlab("Time") 



RE_onepatch_noM_DF_01 = subset(RE_onepatch_DF_noM,
                              RE_onepatch_DF_noM$time > 9100  & 
                                RE_onepatch_DF_noM$time < 9500 & 
                                RE_onepatch_DF_noM$primary_removal == 0.1 )

plot(RE_onepatch_noM_DF_01$H_S,RE_onepatch_noM_DF_01 $N_P,type = 'l')


# Dynamically:
# When you disturb this one-patch system at its equilibrium,
# this leads to a decrease in the total vector abundance, which will
# of course lead to a decrease in the RE (LEFT)! However, if you look at the
# (RIGHT), a decrease in total vector abundance does not actually lead to a
# decrease in RE. It actually leads to a a drastic rise in RE. So I don't think
# it's purely abundance that drives what we are seeing.

### Thoughts: I think this is saying that it's not the total vector abundance
### that is driving the RE! Because the total vector abundance is the same!
### It's just the epidemiolgoical important vector species is switched. In
### the second plot the secondary vector is being repressed and they're
### very good transmitter.

### At the maximum RE for each grouping, what is the total abundance?

Post_Eq_RE_onepatch_DF <- subset(RE_onepatch_DF, RE_onepatch_DF$time > 9124)
Post_Eq_RE_onepatch_DF_Switch <- subset(RE_onepatch_DF_Switch, RE_onepatch_DF_Switch$time > 9124)
Post_Eq_RE_onepatch_DF_noM <- subset(RE_onepatch_DF_noM, RE_onepatch_DF_noM$time > 9124)

Abundance_max_RE_DF <-
  do.call(
    rbind,
    lapply(
      split(
        Post_Eq_RE_onepatch_DF,
        Post_Eq_RE_onepatch_DF$primary_removal
      ),
      function(x) {
        x[which.max(x$RE), ]
      }
    )
  )

Abundance_max_RE_DF_Switch <- do.call(
  rbind,
  lapply(
    split(
      Post_Eq_RE_onepatch_DF_Switch,
      Post_Eq_RE_onepatch_DF_Switch$primary_removal
    ),
    function(x) {
      x[which.max(x$RE), ]
    }
  )
)

Abundance_max_RE_DF_noM <- do.call(
  rbind,
  lapply(
    split(
      Post_Eq_RE_onepatch_DF_noM,
      Post_Eq_RE_onepatch_DF_noM$primary_removal
    ),
    function(x) {
      x[which.max(x$RE), ]
    }
  )
)


Abundance_max_RE_DF$equlibrium <- ifelse(Abundance_max_RE_DF$primary_removal == 1,
  "Y", "N"
)
Abundance_max_RE_DF_Switch$equlibrium <- ifelse(Abundance_max_RE_DF_Switch$primary_removal == 1,
  "Y", "N"
)
Abundance_max_RE_DF_noM$equlibrium <- ifelse(Abundance_max_RE_DF_noM$primary_removal == 1,
  "Y", "N"
)

### This suggests that the total abundance can overshoot the equilibirum

Total_vector_RE_GG <- ggplot(
  Abundance_max_RE_DF,
  aes(
    x = N_P + N_M, y = RE,
    color = as.factor(primary_removal)
  )
) +
  geom_point(size = 2) +
  geom_point(
    data =
      subset(
        Abundance_max_RE_DF,
        Abundance_max_RE_DF$equlibrium == "Y"
      ),
    color = "black", size = 3
  ) +
  scale_color_viridis(discrete = TRUE, name = "") +
  xlab("Total vector abundance") +
  theme_classic()


Total_vector_RE_Switch_GG <- ggplot(
  Abundance_max_RE_DF_Switch,
  aes(
    x = N_P + N_M, y = RE,
    color = as.factor(primary_removal)
  )
) +
  geom_point(size = 2) +
  geom_point(
    data =
      subset(
        Abundance_max_RE_DF,
        Abundance_max_RE_DF$equlibrium == "Y"
      ),
    color = "black", size = 3
  ) +
  scale_color_viridis(discrete = TRUE, name = "") +
  xlab("Total vector abundance") +
  theme_classic()


Total_vector_RE_noM_GG <-
  ggplot(
    Abundance_max_RE_DF_noM,
    aes(
      x = as.numeric(N_P), y = RE,
      color = as.factor(primary_removal)
    )
  ) +
  geom_point(size = 2) +
  geom_point(
    data =
      subset(
        Abundance_max_RE_DF,
        Abundance_max_RE_DF$equlibrium == "Y"
      ),
    color = "black", size = 3
  ) +
  scale_color_viridis(discrete = TRUE, name = "") +
  xlab("Total vector abundance") +
  theme_classic()

Total_vector_RE_GG + Total_vector_RE_Switch_GG + Total_vector_RE_noM_GG


### What is the total maximum abundance regardless of RE?

max_Abundance_RE_DF <- do.call(
  rbind,
  lapply(
    split(
      Post_Eq_RE_onepatch_DF,
      Post_Eq_RE_onepatch_DF$primary_removal
    ),
    function(x) {
      x[which.max(x$N_P + x$N_M), ]
    }
  )
)

ggplot(
  max_Abundance_RE_DF,
  aes(
    x = N_P + N_M, y = RE,
    color = as.factor(primary_removal)
  )
) +
  geom_point(size = 2) +
  geom_point(
    data = subset(
      max_Abundance_RE_DF,
      max_Abundance_RE_DF$equlibrium == "Y"
    ),
    color = "black", size = 3
  ) +
  scale_color_viridis(discrete = TRUE, name = "") +
  xlab("Total vector abundance") +
  theme_classic()


## At what abundance is RE the maximum?
eq_RE_onepatch <- subset(
  RE_onepatch_DF,
  RE_onepatch_DF$time > 9124
)

eq_RE_onepatch_df <- do.call(
  rbind,
  lapply(
    split(
      eq_RE_onepatch,
      eq_RE_onepatch$primary_removal
    ),
    function(x) x[which.max(x$RE), ]
  )
)

ggplot(
  eq_RE_onepatch_df,
  aes(x = N_P, y = N_M, fill = RE)
) +
  geom_point(size = 4, shape = 21, color = "black") +
  xlab("Total primary vectors") +
  ylab("Total secondary vectors") +
  scale_fill_viridis(option = "rocket") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15)
  )

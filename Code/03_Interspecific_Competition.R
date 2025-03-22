# --------------------------------------------------------------------------------
# Research Question:
# How does the change in the primary to secondary coefficient influence RM
# with the change in the transmission efficiency?
# --------------------------------------------------------------------------------
param_standard <- get_parameters("standard")
c_PM_standard <- param_standard["c_PM"] ## Competition effect of p.vector on s.vector
c_MM_standard <- param_standard["c_MM"]
c_MP_standard <- param_standard["c_MP"]
c_PP_standard <- param_standard["c_PP"]

modifier <- seq(.6766667,2,0.01)
Mortality_P <- c(0.01,0.25, 0.5,0.75)


competition_param <-
  data.frame(expand.grid(
    c_PM = c_PM_standard * modifier,
    mortality_P = Mortality_P
  ))


competition_param_list <- vary_parameter_value(
  param_standard, c( "c_PM", "mortality_P"), competition_param
)



RE_COMPETITION <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("c_PM", "mortality_P"),
    vector_value = competition_param
  ) |>
  Calculate_change_baseline(
    competition_param_list,
    competition_param, "No"
  )

#Note: When you analytically solve for when the primary vector competitively 
#excludes the secondary vector it should be when C_PM is multiplied by 1.5.
#There is a slight numerical error when simulating it so I decided to fix it myself


RE_COMPETITION$standardized <- RE_COMPETITION$c_PM/c_PM_standard
RE_COMPETITION[RE_COMPETITION$standardized > 1.49,]$RE <- 0


panel_1 <- 
  ggplot(RE_COMPETITION, aes(
    x = standardized,
    y = (RE), 
    group = as.factor(1-mortality_P)))+
    geom_line() +
  geom_rect(aes(xmin=1.5, xmax=2, ymin=0, ymax=1.25), fill = 'grey', 
            alpha = 0.01)+
  geom_point(data = subset(RE_COMPETITION,RE_COMPETITION$standardized ==1),
             aes(x = standardized, y = RE), size = 1.5) +
  scale_fill_viridis() + 
  geom_vline(xintercept = c(1.5)) + 

  xlab(expression("Multiplier of primary on secondary competition ("*c[PM]*")")) +
  ylab(expression("Increase from baseline " * R[E]^"*")) +
  theme_classic() + 
  theme(legend.position = 'top', 
        axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 15, color = 'black'))
  
panel_2 <- ggplot(
  RE_COMPETITION , aes(
    x = c_PM/c_PM_standard,
    y = (max_NM), color = as.factor(1-mortality_P),
    group = as.factor(1-mortality_P)))+ 
    geom_line()+
  geom_rect(aes(xmin=1.5, xmax=2, ymin=0, ymax=1.25), fill = 'grey', 
            alpha = 1)+
  geom_point(data = subset(RE_COMPETITION,RE_COMPETITION$standardized ==1),
             aes(x = standardized, y = max_NM), size = 1.5, color = 'black') +
    scale_color_grey() + 
  geom_vline(xintercept = c(1.5))+
  xlab(expression("Modifier of primary on secondary competition ("*c[PM]*")")) +
  ylab(expression("Increase from " * N[M]^"*"))+
  theme_classic() + 
   theme(legend.position = 'none',
         axis.text = element_text(size = 14, color = 'black'),
         axis.title = element_text(size = 15, color = 'black'))


panel_1/panel_2





try<- Isocline_Primary (1* c_PM_standard)
try_P <- Phaseplot_Primary(1* c_PM_standard)

try2<- Isocline_Primary (0.75* c_PM_standard)
try2_P <- Phaseplot_Primary(0.75* c_PM_standard)


try3<- Isocline_Primary (1.25* c_PM_standard)
try3_P<- Phaseplot_Primary (1.25* c_PM_standard)



modifier <- c(0.75,1,1.25)
Mortality_P <- c(0.01)


competition_param <-
  data.frame(expand.grid(
    c_PM = c_PM_standard * modifier,
    mortality_P = Mortality_P
  ))


competition_param_list <- vary_parameter_value(
  param_standard, c( "c_PM", "mortality_P"), competition_param
)



RE_COMPETITION_2 <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("c_PM", "mortality_P"),
    vector_value = competition_param
  ) 

testing<- Calculate_Human_Reff_Expanded(
  RE_COMPETITION_2[[2]],competition_param_list[[2]])

eep<- testing[testing$time > competition_param_list[[2]][["disturbance_time"]]-1,]

testingnub<- Calculate_Human_Reff_Expanded(
  RE_COMPETITION_2[[1]],competition_param_list[[1]])
eep1 <- testingnub[testingnub$time > competition_param_list[[1]][["disturbance_time"]]-1,]



testingup<- Calculate_Human_Reff_Expanded(
  RE_COMPETITION_2[[3]],competition_param_list[[3]])
eep2 <- testingup[testingup$time > competition_param_list[[3]][["disturbance_time"]]-1,]





ggplot() + 
  geom_segment(data = try2[[1]], aes(x = x, xend = xend, y = y, yend= yend), color = "#646198") + 
  geom_segment(data = try2[[2]], aes(x = x, xend = xend, y = y, yend= yend), color = "#D65739")+ 
  geom_point(aes(x = 1104.40420,y =  1006.036)) + 
  theme_classic() + 
  scale_color_viridis() + 
  geom_quiver(data = try2_P, aes(x = Var1, y = Var2, u = slope_P, v = slope_M, color = RE),vecsize = 2, alpha = 0.5) +
  coord_cartesian(xlim = c(0,2000), ylim = c(0, 2000)) + 
  geom_path(data = eep1, aes( x= NP, y= NM), size = 1, alpha =1)+
  xlab("Primary vectors") + 
  ylab("Secondary vectors")  +


ggplot() + 
  geom_segment(data = try[[1]], aes(x = x, xend = xend, y = y, yend= yend), color = "#646198", size = 1) + 
  geom_segment(data = try[[2]], aes(x = x, xend = xend, y = y, yend= yend), color = "#D65739", size = 1)+ 
  geom_quiver(data = try_P, aes(x = Var1, y = Var2, u = slope_P, v = slope_M, color = RE),vecsize = 2, alpha = 0.5) +
  theme_classic() + 
  geom_path(data = eep, aes( x= NP, y= NM), size = 1, alpha =1)+
  scale_color_viridis() + 
  geom_point(aes(x = 1106.63082,y =  672.0430)) + 
  coord_cartesian(xlim = c(0,2000), ylim = c(0, 2000)) + 
  xlab("Primary vectors") + 
  ylab("Secondary vectors") +



  
  ggplot() + 
  geom_segment(data = try3[[1]], aes(x = x, xend = xend, y = y, yend= yend), color = "#646198", size = 1) + 
  geom_segment(data = try3[[2]], aes(x = x, xend = xend, y = y, yend= yend), color = "#D65739", size = 1)+ 
  geom_quiver(data = try3_P, aes(x = Var1, y = Var2, u = slope_P, v = slope_M, color = RE),vecsize = 2, alpha = 0.5) +
  theme_classic() + 
  geom_path(data = eep2, aes( x= NP, y= NM), size = 1, alpha =1)+
  scale_color_viridis() + 
  geom_point(aes(x = 1108.86644,y =  336.7003 )) + 
  coord_cartesian(xlim = c(0,2000), ylim = c(0, 2000)) + 
  xlab("Primary vectors") + 
  ylab("Secondary vectors") + 
  
  plot_layout(guides = 'collect')
  

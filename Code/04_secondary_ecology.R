### Role of secondary vector:

# Retrieve "standard" parameters and set the different values of disturbances
param_standard <- get_parameters("standard")
Mortality_P <- seq(0.01,1, 0.01)
modifier <- seq(0, 2, length = 20)


f_P_standard <- param_standard["f_P"] 
theta_P_standard  <- param_standard["theta_P"]
f_M_standard <- param_standard["f_M"] 
theta_M_standard  <- param_standard["theta_M"]


secondary_param <-
  data.frame(expand.grid(
    f_M = f_M_standard * modifier,
    theta_M = theta_M_standard * modifier,
    mortality_P = 0.01
  ))


secondary_param_list <- vary_parameter_value(
  param_standard, c("f_M", "theta_M", "mortality_P"), 
  secondary_param 
)


RE_SECONDARY <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("f_M", "theta_M", "mortality_P"),
    vector_value = secondary_param
  ) |>
  Calculate_change_baseline(
    secondary_param_list,
    secondary_param, "No"
  )


ggplot(RE_SECONDARY ,aes(x = f_M/f_P_standard, 
       y = theta_M/theta_P_standard, fill = RE))+
  geom_tile() + 
  xlab(expression("Modifier of secondary biting rate (" *f[M]/f[P]*")"))+
  ylab(expression("Modifier of secondary transmission probability (" *theta[M]/theta[P]*")"))+
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_viridis(name = expression("Increase from " * R[0])) + 
  theme(legend.position = "top",
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10))

ggsave(here("Main_Figures", "Figure_5.pdf"), width = 6, height = 6, units = "in")

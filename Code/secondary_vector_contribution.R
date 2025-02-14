param_standard <- get_parameters("standard")
Mortality_P <- seq(0.01, 1, length = 10)
modifier <- seq(0, 2, 0.1)

f_M_standard <- param_standard["f_M"] ## Competition effect of s.vector on p.vector
theta_M_standard <- param_standard["theta_M"] ## Competition effect of p.vector on s.vector
f_P_standard <- param_standard["f_P"]
theta_P_standard <- param_standard["theta_P"]


secondary_param <-
  data.frame(expand.grid(
    f_M = f_M_standard * modifier,
    mortality_P = Mortality_P
  ))

secondary_param_list <- vary_parameter_value(param_standard, c("f_M", "mortality_P"), secondary_param)


RE_SM <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("f_M", "mortality_P"),
    vector_value = secondary_param
  ) |>
  Calculate_change_baseline(
    secondary_param_list,
    secondary_param, "No"
  )


heatmap_SM_A <-
  ggplot(RE_SM, aes(x = f_M / f_P_standard, y = 1 - mortality_P, fill = max_mtoH)) +
  geom_tile(size = 0.8) +
  coord_equal() +
  xlab("Transmission efficacy of secondary vector \ncompared to primary vector") +
  ylab("Proportion of primary vectors removed") +
  scale_fill_viridis(option = "mako", name = "Secondary vector contribution to R0") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.position = "top"
  )


heatmap_SM_A <-
  ggplot(RE_SM, aes(x = f_M / f_P_standard, y = 1 - mortality_P, fill = RE)) +
  geom_tile(size = 0.8) +
  coord_equal() +
  xlab("Transmission efficacy of secondary vector \ncompared to primary vector") +
  ylab("Proportion of primary vectors removed") +
  scale_fill_viridis(option = "mako", name = "Secondary vector contribution to R0") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.position = "top"
  )

RE_SM_2 <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("f_M", "mortality_P"),
    vector_value = secondary_param
  ) |>
  Calculate_change_baseline(
    secondary_param_list,
    secondary_param, "Yes"
  )

RE_SM_2$standardized_FM <- RE_SM_2$f_M / f_M_standard
RE_SM_2_subset <- subset(RE_SM_2, RE_SM_2$mortality_P %in% c(0.01) &
  RE_SM_2$standardized_FM %in%
    c(0.5, 1, 1.5, 2))


secondary_line_B <- ggplot(RE_SM_2_subset, aes(
  x = time - 9125,
  y = max_mtoH,
  linetype = as.factor(f_M / f_P_standard)
)) +
  geom_line() +
  xlab("Time since disturbance") +
  ylab("Secondary vector contribution to R0") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black")
  )

heatmap_SM_A + secondary_line_B

ggsave(here("Figures_Process", "GG_secondary_contribution.pdf"),
  units = "in", width = 10, height = 5
)

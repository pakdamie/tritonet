# Vary the competition coefficient of P on M as well as M on M.

param_standard <- get_parameters("standard")
c_PM_standard <- param_standard["c_PM"] ## Competition effect of p.vector on s.vector
c_MM_standard <- param_standard["c_MM"] ## Competition effect of s.vector on s.vector

modifier <- seq(0.01, 2, 0.10)
Mortality_P <- c(0.01, 0.25, 0.5, 0.75)


competition_param <-
  data.frame(expand.grid(
    c_MM = c_MM_standard * modifier,
    c_PM = c_PM_standard * modifier,
    mortality_P = Mortality_P
  ))


competition_param_list <- vary_parameter_value(
  param_standard, c("c_MM", "c_PM", "mortality_P"), competition_param
)


RE_COMPETITION <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("c_MM", "c_PM", "mortality_P"),
    vector_value = competition_param
  ) |>
  Calculate_change_baseline(
    competition_param_list,
    competition_param, "No"
  )


ggplot(RE_COMPETITION, aes(
  x = c_PM / c_PM_standard,
  y = c_MM / c_MM_standard, fill = RE
)) +
  geom_tile(aes(color = RE)) +
  geom_point(aes(x = 1, y= 1), color = 'white', shape = 3) + 
  facet_wrap(~ c(1 - mortality_P)) +
  scale_fill_viridis(name = expression("Change from " * R[0]^"*")) +
  scale_color_viridis(name = expression("Change from " * R[0]^"*")) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  xlab(expression("Modifier of primary on secondary competition ("*c[PM]*")")) +
  ylab(expression("Modifier of secondary on secondary competition ("*c[MM]*")")) +
 
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14)
  )

ggsave(here("Figures", "RE_DF_cmm_cpm_R0.pdf"),
       width = 8, height = 6, units = "in"
)


ggplot(RE_COMPETITION, aes(
  x = c_PM / c_PM_standard,
  y = c_MM / c_MM_standard, fill = max_NM
)) +
  geom_tile(aes(color = max_NM)) +
  geom_point(aes(x = 1, y= 1), color = 'white', shape = 3) + 
  facet_wrap(~ c(1 - mortality_P)) +
  scale_fill_viridis(name = expression("Change from " * N[M]^"*"),option = 'rocket') +
  scale_color_viridis(name = expression("Change from " * N[M]^"*"), option = 'rocket') +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  xlab(expression("Modifier of primary on secondary competition ("*c[PM]*")")) +
  ylab(expression("Modifier of secondary on secondary competition ("*c[MM]*")")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14)
  )

ggsave(here("Figures", "RE_DF_cmm_cpm_NM.pdf"),
  width = 8, height = 6, units = "in"
)

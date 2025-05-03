# VARY THE INTERSPECIFIC P on M as well as M_M
# FIXING MP AS WELL AS PP.
param_standard <- get_parameters("standard")
c_MM_standard <- param_standard["c_MM"] ## Competition effect of p.vector on s.vector

modifier <- seq(0.01, 2, length = 100)
Mortality_P <- c(0.01, 0.25, 0.5, 0.75)

competition_param <-
  data.frame(expand.grid(
    c_MM = c_MM_standard * modifier,
    mortality_P = Mortality_P
  ))


competition_param_list <- vary_parameter_value(
  param_standard, c("c_MM", "mortality_P"), competition_param
)


sub_RE_COMPETITION <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("c_MM", "mortality_P"),
    vector_value = competition_param
  ) |>
  Calculate_change_baseline(
    competition_param_list,
    competition_param, "No"
  )

panel_1 <- ggplot(
  sub_RE_COMPETITION, aes(
    x = c_MM / c_MM_standard,
    y = (RE), color = as.factor(1 - mortality_P),
    group = as.factor(1 - mortality_P)
  )
) +
  geom_line(size = 0.8) +
  scale_color_discrete_sequential(name = "Control\nintensity", palette = "ag_Sunset", rev = FALSE, n = 4) +
  xlab(expression("Multiplier of secondary on secondary competition (" * c[PM] * ")")) +
  ylab(expression("Increase from " * R[E]^"*")) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, color = "black")
  )




panel_1


ggsave(here("Main_Figures", "SUPP_RE_secondarycompetition.pdf"),
  width = 6, height = 6, units = "in"
)

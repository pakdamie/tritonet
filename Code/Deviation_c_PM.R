# VARY THE INTERSPECIFIC P on M as well as M_M
# FIXING MP AS WELL AS PP.
param_standard <- get_parameters("standard")
c_PM_standard <- param_standard["c_PM"] ## Competition effect of p.vector on s.vector

modifier <- seq(0.01, 2, length = 25)
Mortality_P <- seq(0.01,1, 0.01)


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

panel_1 <- ggplot(
  RE_COMPETITION, aes(
    y = 1-mortality_P,
    x = c_PM/c_PM_standard, fill = RE)) + 
  geom_tile() + 
  scale_fill_viridis() + 
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(0,2)) + 
  scale_y_continuous(expand = c(0,0)) + 
  xlab(expression("Modifier of primary on secondary competition ("*c[PM]*")")) +
  ylab("Proportion of primary vectors removed") + 
  theme_classic() + 
  theme(legend.position = 'top', 
        axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 15, color = 'black'))
  
sub_RE_COMPETITION <- subset(RE_COMPETITION,
                             RE_COMPETITION$mortality_P %in% c(0.01,0.25, 0.5,0.75))


panel_2 <- ggplot(
  sub_RE_COMPETITION , aes(
    x = c_PM/c_PM_standard,
    y = (max_NM), color = as.factor(1-mortality_P),
    group = as.factor(1-mortality_P)))+ 
    geom_line(size = 1.2)+
    scale_color_grey() + 
  coord_cartesian(xlim=c(0,2)) + 
  scale_x_continuous(expand = c(0,0))+
  xlab(expression("Modifier of primary on secondary competition ("*c[PM]*")")) +
  ylab(expression("Increase from " * N[M]^"*"))+
  theme_classic() + 
   theme(legend.position = 'none',
         axis.text = element_text(size = 14, color = 'black'),
         axis.title = element_text(size = 15, color = 'black'))


panel_1/panel_2


ggsave(here("Figures", "RE_DF_primarycompetition.pdf"),
  width = 6, height = 10, units = "in"
)


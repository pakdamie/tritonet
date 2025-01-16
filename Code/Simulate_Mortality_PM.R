### Supplementary code/figure

# How does the disturbance intensity of the primary and secondary vectors influence RE?


# Retrieve "standard" parameters and set the different values of disturbances
param_standard <- get_parameters("standard")
Mort_P <- seq(0.01, 1, 0.1) # primary
Mort_M <- seq(0.01, 1, 0.1) # secondary

Mort_PM <- expand.grid("mortality_P" = Mort_P, "mortality_M" = Mort_M)

# Simulate model output and calculate the RE
RE_mortality_PM <-
  Simulate_Model_Output(
    param_standard,
    c("mortality_P", "mortality_M"), Mort_PM
  ) |>
  lapply(Calculate_Human_REff, param = param_standard)

# Assign mortality_P and mortality_M value to each list element
for (i in seq_along(RE_mortality_PM)) {
  RE_mortality_PM[[i]]$mort_P <- Mort_PM[i, "mortality_P"]
  RE_mortality_PM[[i]]$mort_M <- Mort_PM[i, "mortality_M"]
}
RE_mortality_PM_DF <- do.call(rbind, RE_mortality_PM)

### Subset time period of interest (slightly before disturbance and post disturbance)
RE_mortality_PM_post <- subset(
  RE_mortality_PM_DF,
  RE_mortality_PM_DF$time > 8900 &
    RE_mortality_PM_DF$time < 9125 + 5000
)

maximum_PM <- aggregate(RE ~ mort_P + mort_M,
  data = RE_mortality_PM_post,
  FUN = function(x) max(x)
)

ggplot(
  maximum_PM,
  aes(x = 1 - mort_P, y = 1 - mort_M, fill = RE)
) +
  geom_tile() +
  xlab("Proportion of primary vector removed") +
  ylab("Proportion of secondary vector removed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_continuous_divergingx(
    name = expression(R[E]),
    mid = 1, n_interp = 21, palette = "Roma", rev = TRUE
  ) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15)
  ) +
  coord_equal()

ggsave(here("Figures", "supp_pm_removal.pdf"),
  units = "in", width = 11, height = 11
)

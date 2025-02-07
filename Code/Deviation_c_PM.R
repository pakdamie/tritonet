# VARY THE INTERSPECIFIC P on M as well as M_M
# But fix the MP AS WELL AS PP.
# Fixed:
c_MP_standard <- 3e-6 ## Competition effect of s.vector on p.vector
c_PP_standard <- 4.5e-4 ## Competition effect of p.vector on s.vector
c_PM_standard <- 3e-4 ## Competition effect of p.vector on s.vector
c_MM_standard <- 2.5e-4 ## Competition effect of s.vector on s.vector

modifier <- seq(0.01, 2, 0.01)
Mortality_P <- c(0.01, 0.25, 0.5, 0.75)


competition_param <-
  data.frame(expand.grid(
    c_PP = c_PP_standard,
    c_MM = c_MM_standard,
    c_MP = c_MP_standard,
    c_PM = c_PM_standard * modifier,
    mortality_P = Mortality_P
  ))


# stupid way that i'm gonna fix

param_standard_list <- NULL
for (i in seq(1:nrow(competition_param))) {
  param_copy <- get_parameters("standard")
  param_copy["c_PP"] <- competition_param[i, "c_PP"]
  param_copy["c_MM"] <- competition_param[i, "c_MM"]
  param_copy["c_PM"] <- competition_param[i, "c_PM"]
  param_copy["c_MP"] <- competition_param[i, "c_MP"]
  param_copy["mortality_P"] <- competition_param[i, "mortality_P"]
  param_standard_list[[i]] <- param_copy
}


RE_CM <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("c_PP", "C_MM", "c_MP", "c_PM", "mortality_P"),
    vector_value = competition_param
  )


RE_DF_inter <- NULL
for (k in 1:length(RE_CM)) {
  RE_tmp <- Calculate_Human_Reff_Expanded(RE_CM[[k]], param_standard_list[[k]])

  eq_RE <- RE_tmp[9124, ]$RE
  eq_NM <- RE_tmp[9124, ]$NM
  eq_NP <- RE_tmp[9124, ]$NP

  post_df <- na.omit(subset(RE_tmp, RE_tmp$time > 9124))

  RE_max <- cbind(competition_param[k, ],
    max_NP = max(post_df$NP) - eq_NP,
    max_NM = max(post_df$NM) - eq_NM,
    min_NP = max(post_df$NP) - eq_NP,
    min_NM = max(post_df$NM),
    max_NV = max(post_df$NM + post_df$NP) - (eq_NM + eq_NP),
    RE = max(post_df$RE) - eq_RE,
    min_RE = min(post_df$RE)
  )

  RE_DF_inter[[k]] <- RE_max
}

RE_DF_inter <- do.call(rbind, RE_DF_inter)

RE_DF_inter$ext_NP <- ifelse(RE_DF_inter$min_NP < 1, "Extinct", "Not")
RE_DF_inter$ext_NM <- ifelse(RE_DF_inter$min_NM < 1, "Extinct", "Not")
RE_DF_inter$coexistence <- ifelse(RE_DF_inter$ext_NP ==
  RE_DF_inter$ext_NM, "coexistence", "Not")


ggplot(RE_DF_inter, aes(
  x = c_PM / c_MP_standard, y = RE, color = max_NM,
  group = as.factor(1 - mortality_P)
)) +
  geom_path(size = 2, linejoin = "mitre", lineend = "round") +
  xlab("Modifier of primary competition") +
  ylab(expression(paste("Deviation from the", R[0]^"*"))) +
  scale_color_viridis(
    name = expression(paste("Deviation from\nthe", N[P]^"*")),
    option = "rocket"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black")
  )


ggsave(here("Figures", "RE_DF_primarycompetition.pdf"),
  width = 8, height = 6, units = "in"
)

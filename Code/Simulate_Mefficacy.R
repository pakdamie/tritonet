###NEED TO BE CLEANED UP ALOT.


### Calculate recovery time
Mort_P <- c(0.01, 0.25, 0.5, 0.75)
lambda_modifier <- seq(0.75, 1.25, 0.01)
f_M_standard <- get_parameters("standard")["f_M"]
theta_M_standard <- get_parameters("standard")["theta_M"]
lambda_modified_grid <- expand.grid(lambda_modifier, Mort_P)
f_M_vec <- f_M_standard * lambda_modified_grid[, 1]
theta_M_vec <- theta_M_standard * lambda_modified_grid[, 1]

secondary_param <- cbind.data.frame(
  f_M = f_M_vec,
  theta_M = theta_M_vec,
  mortality_P = lambda_modified_grid[, 2]
)

simulated_model <- Simulate_Model_Output_PostD(c("f_M", "theta_M", "mortality_P"), secondary_param)

param_standard_list <- NULL
for (i in seq(1:nrow(secondary_param))) {
  param_copy <- get_parameters("post_disturb")

  param_copy["f_M"] <- secondary_param[i, "f_M"]
  param_copy["theta_M"] <- secondary_param[i, "theta_M"]
  param_copy["mortality_P"] <- secondary_param[i, "mortality_P"]
  param_standard_list[[i]] <- param_copy
}


RE_DF_waitingtime <- NULL
RE_DF_maxRE <- NULL

for (k in 1:length(simulated_model)) {
  RE_tmp <- Calculate_Human_REff(
    simulated_model[[k]],
    param_standard_list[[k]]
  )


  eq <- 0.9579212
  full_df <- cbind.data.frame(time = c(0, RE_tmp$time), RE = c(eq, RE_tmp$RE))
  steady_state <- find_steady_state(full_df$RE, max_diff = 10^(-5))[1]


  RE_DF_waitingtime[[k]] <- steady_state
  RE_DF_maxRE[[k]] <- full_df[which.max(full_df$RE), ]$time
}


RE_DF_waitingtime <- cbind(lambda_modified_grid,
  recovery_time = do.call(rbind, RE_DF_waitingtime)
)

RE_DF_maxRE <- cbind(lambda_modified_grid,
  maxRE = do.call(rbind, RE_DF_maxRE)
)

waiting_time_GG <- ggplot(RE_DF_waitingtime, aes(
  x = Var1, y = recovery_time,
  color = as.factor(1 - Var2),
  group = as.factor(1 - Var2)
)) +
  geom_line(size = 1.2) +
  xlab(expression("Multiplier of " * lambda[M])) +
  ylab("Time to recovery") +
  scale_color_viridis(discrete = TRUE, name = "Disturbance severity") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black")
  )


max_RE_GG <- ggplot(RE_DF_maxRE, aes(
  x = Var1, y = maxRE,
  color = as.factor(1 - Var2),
  group = as.factor(1 - Var2)
)) +
  geom_line(size = 1.2) +
  xlab(expression("Multiplier of " * lambda[M])) +
  ylab(expression("Time to maximum " * R[E])) +
  scale_color_viridis(discrete = TRUE, name = "Disturbance severity") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black")
  )

max_RE_GG / waiting_time_GG + plot_layout(guides = "collect")

### Role of secondary vector:

# Retrieve "standard" parameters and set the different values of disturbances
param_standard <- get_parameters("standard")
Mortality_P <- data.frame(Mort_P = c(0.25))
modifier <- seq(0, 2, length = 20)

f_M_standard <- 0.25 * 0.75 ## Competition effect of s.vector on p.vector
theta_M_standard <- 0.50 * 0.75 ## Competition effect of p.vector on s.vector

f_P_standard <- 0.25
theta_P_standard <-  0.50


secondary_param <-
  data.frame(expand.grid(
    f_M = f_M_standard * modifier,
    theta_M = theta_M_standard * modifier,
    mortality_P = 0
  ))


# stupid way that i'm gonna fix

param_standard_list <- NULL
for (i in seq(1:nrow(secondary_param))) {
        param_copy <- get_parameters("standard")
        param_copy["f_M"] <- secondary_param[i, "f_M"]
        param_copy["theta_M"] <- secondary_param[i, "theta_M"]
        param_standard_list[[i]] <- param_copy
}


RE_SM <- Simulate_Model_Output(
                parameter = get_parameters("standard"),
                infection_start = "No",
                variable_interest = c("f_M","theta_M", "mortality_P"),
                vector_value = secondary_param
        )



RE_DF_second <- NULL
for (k in 1:length(RE_SM)) {
        RE_tmp <- Calculate_Human_REff(RE_SM[[k]], param_standard_list[[k]])
        eq_RE <- RE_tmp[9124, ]$RE
        eq_NM <- RE_tmp[9124, ]$NM
        eq_NP <- RE_tmp[9124, ]$NP
        
        post_df <- na.omit(subset(RE_tmp, RE_tmp$time > 9124))
        
        RE_max <- cbind(secondary_param [k, ],
                        max_NP = max(post_df$NP) - eq_NP,
                        max_NM = max(post_df$NM) - eq_NM,
                        min_NP = max(post_df$NP) - eq_NP,
                        min_NM = max(post_df$NM),
                        max_NV = max(post_df$NM + post_df$NP) - (eq_NM + eq_NP),
                        RE = max(post_df$RE) -eq_RE ,
                        min_RE = min(post_df$RE)
        )
        
        RE_DF_second[[k]] <- RE_max
}

RE_DF_seconddf <- do.call(rbind, RE_DF_second )

ggplot(RE_DF_seconddf, aes(x =  f_M/f_P_standard, y= theta_M/theta_P_standard, fill = RE)) +
        geom_tile() + scale_fill_viridis(option = "rocket") + 
        


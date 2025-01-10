f_m_standard <- get_parameters("no_diff")["f_M"]
theta_m_standard <- get_parameters("no_diff")["theta_M"]


# Interspecific competition
modifier_values <- seq(0.05, 2.5, 0.10)
modified_theta_M <- as.numeric(modifier_values * get_parameters("no_diff")[["theta_M"]])
modified_f_M <- as.numeric(modifier_values * get_parameters("no_diff")[["f_M"]])
Mortality_P <- 0.01
thetaf_M_disturbance <- expand.grid(theta_M = modified_theta_M, 
                                    f_M = modified_f_M,
                                    mortality_P = Mortality_P)

# stupid way that i'm gonna fix

param_standard_list <- NULL
for (i in seq(1:nrow(thetaf_M_disturbance))) {
        param_copy <- get_parameters("standard")
        
        param_copy["theta_M "] <- thetaf_M_disturbance[i, "theta_M"]
        param_copy["f_M"] <- thetaf_M_disturbance[i, "f_M"]
        param_copy["mortality_P"] <- thetaf_M_disturbance[i, "mortality_P"]
        param_standard_list[[i]] <- param_copy
}


RE_CM <-
        Simulate_Model_Output(
                get_parameters("standard"),
                c("theta_M", "f_M", "mortality_P"),
                thetaf_M_disturbance
        )



RE_DF_fm_theta <- NULL
for (k in 1:length(RE_CM)) {
        RE_tmp <- Calculate_Human_REff(RE_CM[[k]], param_standard_list[[k]])
        post_df <- na.omit(subset(RE_tmp, RE_tmp$time > 9124))
        
        if (nrow(post_df) == 0) {
                RE_max <- cbind(thetaf_M_disturbance[k, ], RE = NA, MtoH = NA)
        } else {
                RE_max <- cbind(thetaf_M_disturbance[k, ], RE = max(post_df$RE),  
                                MtoH = post_df[which.max(post_df$MtoH),]$MtoH/
                                        post_df[which.max(post_df$MtoH),]$RE)
        }
        RE_DF_fm_theta [[k]] <- RE_max
}

RE_DF_fm_theta  <- do.call(rbind, RE_DF_fm_theta )

 ggplot(
        RE_DF_fm_theta,
        aes(
                x = as.factor(theta_M/theta_m_standard),
                y = as.factor(f_M/f_m_standard),
                fill = RE
        )
) +
        geom_tile() +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        ggtitle("A. Interspecific competition") + 
         scale_fill_binned_divergingx(
                name = expression(R[E]),
                mid = 1, n_interp = 11,
                palette = "Roma", rev = TRUE, breaks = c(0.5,1,5,10,20)
        ) +
        xlab("Transmission probability") +
        ylab("Biting rate")+ 
        theme(axis.text = element_text(size = 14, color = 'black'),
              axis.title = element_text(size = 14, color = 'black')) 



 RE_DF_fm_theta$M <- RE_DF_fm_theta$theta_M * RE_DF_fm_theta$f_M
 RE_DF_fm_theta$Above1 <- ifelse( RE_DF_fm_theta$RE >= 1, "True", "False")
 

 
mtoH_inter_GG <- ggplot(
        RE_DF_inter,
        aes(
                x = as.factor(c_PM / c_PM_standard),
                y = as.factor(c_MP / c_MP_standard),
                fill = MtoH
        )
) +
        geom_tile() +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        ggtitle("A. Interspecific competition") + 
        scale_fill_viridis() +
        xlab(expression("Modifier of primary on secondary competition"~(italic(c[PM])))) +
        ylab(expression("Modifier of secondary on primary competition"~(italic(c[MP])))) + 
        theme(axis.text = element_text(size = 14, color = 'black'),
              axis.title = element_text(size = 14, color = 'black')) 


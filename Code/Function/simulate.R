###
mid_connectance_network <- simulate_spatial_network(24601, 20, 0.10)

impact_on_primary <- seq(0,1,0.1)
spatial_coverage <- 0.5

expanded_param_test <- expand.grid(
        impact_P = impact_on_primary, 
        sp_cov = spatial_coverage)




###The function to simulate a lot of infection
simulate_infection <- function(params) {
    impact_P <- params$impact_P
    sp_cov <- params$sp_cov

return(discrete_trito_model(1000, param_ALT, 500, impact_P, 0.9, sp_cov, 
                                             mid_connectance_network , 1)
                )
        }
        
results <- apply(expanded_param_test, 1, function(row) {
    params <- list(
        impact_P = as.numeric(row["impact_P"]), 
        sp_cov = as.numeric(row["sp_cov"])
    )
    simulate_infection(params)
})

RE <- do.call(rbind,unlist(lapply(results, function(x) calculate_R_effective_discrete_patch(param_ALT,x, disturbance = 500)),
             recursive = FALSE))

max_Primary <-do.call(rbind, lapply(results, function(x) {
        sum_row <- sum(x[[4]][501, ] + x[[5]][501, ] + x[[6]][501, ] + x[[7]][501, ])
        primary <- sum(x[[4]][501, ] + x[[5]][501, ])
        secondary <- sum(x[[6]][501, ] + x[[7]][501, ])
        return(data.frame(sum_row, primary,secondary))
}))

 RE <- cbind(RE,max_Primary)

RE$ratio <- RE$primary/((RE$primary + RE$secondary))
RE <- cbind(RE, expanded_param_test )

ggplot(RE, aes(x = 1-impact_P , y = RE, color = sum_row)) + 
    geom_point(size = 8) +
    xlab("Proportion of p.species removed") + 
    ylab("Maximum effective reproductive number") + 
    theme_classic() + 
    theme(axis.text = element_text(size = 14, color = 'black'),
          axis.title = element_text(size = 15, color = 'black')) + 
    scale_color_viridis()

ggplot(RE, aes(x = 1-impact_P , y = RE, color = ratio)) + 
    geom_point(size = 8) +
    xlab("Proportion of p.species removed") + 
    ylab("Maximum effective reproductive number") + 
    theme_classic() + 
    theme(axis.text = element_text(size = 14, color = 'black'),
          axis.title = element_text(size = 15, color = 'black')) + 
    scale_color_viridis()


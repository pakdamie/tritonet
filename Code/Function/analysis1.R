param_standard <-
  c(b_H = 1 / 1000, ## Human mortality rate
    b_P = 0.01, # P. Vector birth rate
    b_M = 0.01, # S. Vector birth rate
    mu_H = 1 / 1000, ## Human death rate
    f_P = 0.02, # Biting rate of the p. vector
    f_M = 0.02 * 0.75, # Biting rate of the s.vector
    theta_P = 0.7, # Transmission probability of p. vector
    theta_M = 0.7 * 0.75, # Transmission probability of s. vector
    theta_H = 0.5, # Transmission probability of human
    gamma = 1 / 90, # Recovery rate of infected human
    c_PM = 4e-6, ## Competition effect of p.vector on s.vector
    c_MP = 2e-6, ## Competition effect of s.vector on p.vector
    c_PP = 4.5e-6, ## Competition effect of p.vector on s.vector
    c_MM = 3e-6, ## Competition effect of s.vector on s.vector
    ntime = 365 * 50, # How long to run the simulation for
    disturbance_time = 365 * 25, # When to disturb the system
    delta_T = 1, # Time-step 1
    mortality_P = 0.01, # This will change
    mortality_M = 1
  )

Initial_List <- create_initial_states(param_standard, patch_num = 1)

model_output_mort_P <- discrete_trito_model_rcpp_ONEPATCH(
        HS = Initial_List[[1]],
        HI = Initial_List[[2]],
        HR = Initial_List[[3]],
        PS = Initial_List[[4]],
        PI = Initial_List[[5]],
        MS = Initial_List[[6]],
        MI = Initial_List[[7]],
        param = param_standard
)


RE_ANALYTICAL <-
        Calculate_analytical_REff(
        model_output_mort_P, # Model outputs
        param_standard
)










mortality_P <- seq(0.1, 1, 0.10)

create_param_list <- function(mortality_value, param_base) {
        param_base_copy <- param_base
        param_base_copy["mortality_P"] <- mortality_value
        return(param_base_copy)

standard_param_L <- lapply(mortality_P, function(m) create_param_list(m, param_standard))



## Damie: Need to rewrite function to include LIST instead of subsetting
model_output_mort_P <- lapply(standard_param_L, function(param) {
        discrete_trito_model_rcpp_ONEPATCH(
                HS = Initial_List[[1]],
                HI = Initial_List[[2]],
                HR = Initial_List[[3]],
                PS = Initial_List[[4]],
                PI = Initial_List[[5]],
                MS = Initial_List[[6]],
                MI = Initial_List[[7]],
                param = param
        )
})

RE_onepatch_List <- NULL
## figure out how to simplify this
for (p in seq(1, length(standard_param_L))) {
        RE_onepatch_List[[p]] <- cbind.data.frame(
                calculate_RE_onepatch(
                        model_output_mort_P[[p]], # Model outputs
                        standard_param_L[[p]], 1
                ),
                primary_removal = standard_param_L[[p]]["mortality_P"]
        )
}

RE_onepatch_DF <- do.call(rbind, RE_onepatch_List)
colnames(RE_onepatch_DF) <- c(
        "RE", "N_P", "N_M", "H_S", "P_S", "M_S",
        "H_I", "P_I", "M_I",
        "time", "primary_removal"
)

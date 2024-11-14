param_standard <-  
        c(b_H = 1/(1000), ##Human mortality rate
          b_P = 0.01, # P. Vector birth rate
          b_M = 0.01,  # S. Vector birth rate
          mu_H =1/(1000),  ##Human death rate
          mu_V = 0,  ##P. vector death rate
          f_P = 0.050, # Biting rate of the p. vector
          f_M = 0.040, # Biting rate of the s.vector
          theta_P = 0.90, # Transmission probability of p. vector
          theta_M = 0.90, # Transmission probability of s. vector
          theta_H  = 0.50, # Transmission probability of human
          gamma = 1/100,  ##Recovery rate of infected human
          c_PM = 5e-6, ##Competition effect of p.vector on s.vector
          c_MP = 1e-6,  ##Competition effect of s.vector on p.vector
          c_PP = 7e-6, ##Competition effect of p.vector on s.vector
          c_MM = 4e-6,
          d = 0.10,
          ntime = 10000,
          disturbance_time = 5000,
          delta_T = 1,
          prop = 1,
          mortality_P = 0.10,
          mortality_M = 0.95) # Transmission probability of s. vector)




adjacency_matrix <- as_adjacency_matrix(
        low_connectance_network,
        type = "both",
        attr = "weight",
        names = TRUE,
        sparse = FALSE
)

two_patch_adjacency_matrix <- matrix(c(0,1,1,0), nrow = 2, ncol = 2,byrow = TRUE)


summed_prob <- rowSums(adjacency_matrix)
adjacency_matrix_adj <- sweep(adjacency_matrix, MARGIN = 1, summed_prob, `/`)

patch_num <- nrow(two_patch_adjacency_matrix)


compartment_label <- c(
        "HS_mat", "HI_mat", "HR_mat",  # Humans (sus/inf/rec)
        "PS_mat", "PI_mat",  # Primary (sus/inf)
        "SS_mat", "SI_mat"   # Secondary(sus/inf)
)

# Using the above compartment label, we then create an empty matrix where 
# the row are the time-steps and the columns are the individual patches

for (i in 1:length(compartment_label)) {
        assign(
                compartment_label[i],
                matrix(0, nrow = param_ALT['ntime'], ncol = patch_num)
        )
}

#Initial conditions (everyone starts out with the same number)
HS_mat[1, ] <- rep(1000, patch_num) 
HI_mat[1, ] <- rep(0, patch_num) 
HR_mat[1, ] <- rep(0, patch_num)

PS_mat[1, ] <- rep(1000, patch_num)
PI_mat[1, ] <- rep(10, patch_num)
SS_mat[1, ] <- rep(1000, patch_num)
SI_mat[1, ] <- rep(10, patch_num)


test <- discrete_trito_model_rcpp(HS = HS_mat,
                          HI = HI_mat,
                          HR = HR_mat, 
                          PS = PS_mat,
                          PI = PI_mat,
                          SS = SS_mat,
                          SI = SI_mat,
                          adj = adjacency_matrix_adj,
                          param_ALT)    

plot_list_groups(test[1:7])[[1]]

tmp <- calculate_R_effective_discrete_patch(param_ALT, test, 5000)


ggplot(tmp[[1]], aes( x= time, 
                      y= RE, 
                      color = patch_num, group = patch_num)) + 
        geom_line() + 
        scale_color_viridis()


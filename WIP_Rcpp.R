low_connectance_network<- simulate_spatial_network(24601, 20, 0.10)

param_standard <-  
  c(b_H = 1/(1000), ##Human mortality rate
    b_P = 0.1, # P. Vector birth rate
    b_M = 0.1,  # S. Vector birth rate
    mu_H =1/(1000),  ##Human death rate
    f_P = 0.05, # Biting rate of the p. vector
    f_M = 0.025, # Biting rate of the s.vector
    theta_P = 0.90, # Transmission probability of p. vector
    theta_M = 0.75, # Transmission probability of s. vector
    theta_H  = 0.50, # Transmission probability of human
    gamma = 1/100 ,  # Recovery rate of infected human
    c_PM = 2e-5, ## Competition effect of p.vector on s.vector
    c_MP = 1e-5, ## Competition effect of s.vector on p.vector
    c_PP = 3e-5, ## Competition effect of p.vector on s.vector
    c_MM = 3e-5, ## Competition effect of s.vector on s.vector
    ntime = 365 * 50,
    disturbance_time = 365 * 25 ,
    delta_T = 1,
    d = 0.05,
    prop = 0.25,
    mortality_P = 0.25,
    mortality_M = 0.95) 

adjacency_matrix <- as_adjacency_matrix(
        low_connectance_network,
        type = "both",
        attr = "weight",
        names = TRUE,
        sparse = FALSE
)

summed_prob <- rowSums(adjacency_matrix)
adjacency_matrix_adj <- sweep(adjacency_matrix, MARGIN = 1, summed_prob, `/`)
patch_num <- nrow(adjacency_matrix_adj )
initial_states<- create_initial_states(param_standard,patch_num)


model<- discrete_trito_model_rcpp(
        HS = initial_states[[1]],
        HI = initial_states[[2]],
        HR = initial_states[[3]], 
        PS = initial_states[[4]],
        PI = initial_states[[5]],
        MS = initial_states[[6]],
        MI = initial_states[[7]],
        adj = adjacency_matrix_adj,
        param_standard)    



plot_list_groups(model[1:7])

calculated_R_patch <- 
  calculate_R_effective_discrete_patch(param_standard, model, 8000 )[[1]]





calculated_R_net <- calculate_R_effective_discrete_net(param_standard, model, 8000)





patch_interest <- model[[8]] + 1
calculated_R_patch$targeted = NA

calculated_R_patch$targeted<- ifelse(calculated_R_patch$patch_num 
                                  %in% patch_interest, "Yes", "No")

calculated_R_patch_2 <- do.call(rbind,lapply(split(calculated_R_patch, calculated_R_patch$patch_num),
       function(x) {x$diff <- x$RE - x$RE[1] 
       return(x)}))


RE_panel_A <- 
ggplot(calculated_R_patch,
       aes( x= time, y= RE, 
            color = targeted, 
            group = patch_num)) + 
        geom_line(alpha = 0.95, size = 0.8) + 
        scale_color_manual(name = "Target",
                           values = c("No" = "darkgrey", "Yes" = "red")) + 
        xlab("Time") + 
        ylab("RE") + 
        geom_vline(xintercept = 8000)+
        theme_classic() + 
        theme(axis.text = element_text(size = 14), 
              axis.title = element_text(size = 15))


RE_panel_B <- 
   ggplot(calculated_R_patch_2,
       aes( x= time, y= diff, 
            color = targeted, 
            group = patch_num)) + 
        geom_line(alpha = 0.95, size = 0.8) + 
        scale_color_manual(name = "Target",values = c("No" = "darkgrey", 
                                      "Yes" = "red")) + 
        xlab("Time") + 
        ylab("Deviation from equilibirum RE") + 
        geom_vline(xintercept = 8000)+
        theme_classic() + 
        theme(axis.text = element_text(size = 14), 
              axis.title = element_text(size = 15))


ggplot(calculated_R_net, 
       aes(x = time, y = RE)) + geom_line(()



RE_panel_A /RE_panel_B + plot_layout(guides = 'collect')

patch_interest <- model[[8]] + 1
nodes_list <- seq(1,100)
V(low_connectance_network)$target <-  ifelse(nodes_list  
                                     %in% patch_interest, "Yes", "No")


new_g2 <- ggnetwork(low_connectance_network, layout = igraph::layout.reingold.tilford(
        low_connectance_network, circular = TRUE))
ggplot(new_g2, 
       aes(x = x, y = y, 
           xend = xend, yend = yend)) +
        geom_edges(size = 0.5, color = 'grey') +
        geom_nodes(aes(color = target), size = 5) + 
        geom_nodetext(aes(label = name), size = 2) +
        scale_color_manual(values = c("No" = "grey", "Yes" = "pink")) + 
        theme_blank()

### 49 Is a strange one
degree(low_connectance_network)[49] #15 edges
neighbors_49 <- neighbors(low_connectance_network, 49, mode = "all")


calculated_R_patch <- 
        calculate_R_effective_discrete_patch(param_standard, model, 8000 )[[1]]


calculated_R_patch$targeted<- ifelse(calculated_R_patch$patch_num 
                                     %in% patch_interest, "Yes", "No")

calculated_R_patch_2 <- do.call(rbind,lapply(split(calculated_R_patch, calculated_R_patch$patch_num),
                                             function(x) {x$diff <- x$RE - x$RE[1] 
                                             return(x)}))

neighbors_49_RE <- subset(calculated_R_patch_2,
                          calculated_R_patch_2 $patch_num %in%
                                  c(49, neighbors_49))

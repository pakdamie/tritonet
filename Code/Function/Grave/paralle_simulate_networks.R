num_cores <- detectCores() - 1  # Use all cores minus one
connectance_interest <- seq(0.05, 0.5, 0.05)

igraph_list <- NULL

simulate_networks <- function(connectance_value) {
        seed_no <- sample(seq(1, 24601, 1), 1)
        igraph_list_sub <- replicate(
                10,
                simulate_spatial_network(seed_no, 20, connectance_value)
        )
        return(igraph_list_sub)
}

igraph_list <- mclapply(connectance_interest, simulate_networks, mc.cores = num_cores)


save(igraph_list , file = "igraph_list.RData")



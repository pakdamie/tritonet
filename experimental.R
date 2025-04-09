###How does connectance, coverage, and intensity influence the maximum RE of the total network


connectance_networks <- replicate(5, simulate_spatial_network(30,5,c(0.15,0.20,0.30)))


replicated_networks <- function(rep_networks){
  
  col_list = NULL
  for (i in 1:ncol(connectance_networks)){
    
    rep_interest <- connectance_networks[,i]
    graph_objects <- rep_interest[[1]]
    
    #Number of nodes 
    node_number <- unique(rep_interest[[2]]$num_nodes)
    
    Initial_list <- create_initial_states(standard_param,
                                          patch_num = node_number)
    
    adjacency_matrices <- lapply(graph_objects,
                                 function(x) 
                                   simulate_as_adjacency_matrix(x))
    
    col_list[[i]] = list(Initial_list, adjacency_matrices)
  }
  return(col_list)
}



rep_list <- replicated_networks(connectance_networks)

sim_model <- function(networks, parameter_list){
  
  
  model_output_big_list <- NULL
  
  for (i in seq(length(networks))){
    
    Initial_list<- rep_list[[i]][[1]]
    Adj_matrix <- rep_list[[i]][[2]]
    
    model_output_list <- NULL
    
    for (p in seq(1, length(parameter_list))) {
      model_output_list[[p]] <-
        two_vector_one_host_disp(
          HS = Initial_list[[1]],
          HI = Initial_list[[2]],
          HR = Initial_list[[3]],
          PS = Initial_list[[4]],
          PI = Initial_list[[5]],
          MS = Initial_list[[6]],
          MI = Initial_list[[7]],
          adj = Adj_matrix[[i]],
          param = parameter_list[[p]]
        )
    }
    
    model_output_big_list[[i]] =  model_output_list[[p]]
  }
}

sim_model(rep_list,parameter_mortality_P_list)


### Setting initial conditions
Initial_List <- create_initial_states(standard_param, 
                                      patch_num)

low_connectance_network <- simulate_spatial_network(24601, 30, 0.10)

adjacency_matrix <- as_adjacency_matrix(
  low_connectance_network,
  type = "both",
  attr = "weight",
  names = TRUE,
  sparse = FALSE
)
summed_prob <- rowSums(adjacency_matrix)
adjacency_matrix_adj <- sweep(adjacency_matrix, MARGIN = 1, summed_prob, `/`)


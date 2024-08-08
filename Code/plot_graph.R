plot_the_graph

adjacency_matrix <- simulate_final_adjacency_matrix(24601, num_patch,
                                                    connectance,max_distance)

g9 <- graph_from_adjacency_matrix(adjacency_matrix , weighted=TRUE,
                                  mode="plus", diag=FALSE)

#' Simulate the adjacency matrix for the main simulation
#' 
#' As the main simulation requires a spatial network, this
#' function generates the necessary adjacency matrix for
#' dispersal.
#'
#' @param num_patch an integer, number of patches that you want to simulate
#' @param connectance a numeric variable (0-1), the connectance of the network
#' @param seed a numeric variable, the random seed generator (default is 24601)
#' @return a matrix (the adjacency matrix) 
#'
#' @examples simulate_adjacency_matrix(100, 10, 0.2)
#' 

simulate_adjacency_matrix <- function(num_patch, 
                                      max_distance, 
                                      connectance, 
                                      seed = 24601) {
  
 set.seed(seed) #Sets the number of seeds
  
  ###Generate the patches and their coordinates
  xy <- seq(1, max_distance, length.out = 2000) #possible coordinates
  x_coord <- sample(xy, num_patch, replace = TRUE) #x-coordinate
  y_coord <- sample(xy, num_patch, replace = TRUE) #y-coordinate 
  xy_coord <- cbind(x_coord, y_coord) #xy-coordinates 
  NegExpDist <- as.matrix(exp(-dist(xy_coord))) #distance matrice with kernel
  
  ###Create an graph
  Adj_graph <- graph.adjacency(NegExpDist, mode = 'undirected', 
                               diag = FALSE, 
                               weighted = TRUE)
  
  ##Adding latitude and longitude
  V(Adj_graph)$Long <- xy_coord [,1]
  V(Adj_graph)$Lat <-  xy_coord [,2]

  
  ### The number of edges we need to calculate the connectance.
  number_of_edges <- (connectance * (num_patch^2))
  

  ###If the number of edges required for connectance is 10, then 
  ### choose the 10 likeliest (highest weight) edges.
  deleted_edges_graph <- delete.edges(Adj_graph,
                         which(E(Adj_graph)$weight < sort(E(Adj_graph)$weight, 
                         decreasing = T)[number_of_edges])) 
 
  ###This will likely lead to partitions into different components

  ###The graph is split into components if the length of the decomposition is not 1!
  is_it_split <- ifelse(length(decompose(deleted_edges_graph)) != 1, "split", "complete")
  

  while(is_it_split == "split"){
         
         for (c in seq(1, length(component_graph))){
                
             component_graphs <- decompose(deleted_edges_graph)
                 
             vertex_component <-  as_ids(V(component_graphs[[c]])) ###identify vertex in component of interest
             
             ###identify the vertices that are not in the same component
             non_component_vertex = as_ids(all_vertex[!(all_vertex %in%  vertex_component)]) 
             
             ###create possible edges
             possible_edges <- expand.grid(as.numeric(vertex_component), as.numeric(non_component_vertex))
             
             ###Get the weight from the original adjacency matrix 
             possible_edges$weight <- E(Adj_graph)[possible_edges[,1] %--% possible_edges[,2]]$weight
             
             ###Get the maximum weight edge
             maximum_edgeweight <- possible_edges[which.max(  possible_edges$weight),]
             
             ###Delete the lowest-weight edge
             deleted_edge <-  delete.edges(deleted_edges_graph,
                                 which.min(E(deleted_edges_graph)$weight))
             
             ###Add the maximum weight edge 
             deleted_edges_graph <- add_edges(deleted_edge,
                                      c(maximum_edgeweight[,1],
                                        maximum_edgeweight[,2]))
             
             
            ###TAD's ADVICE: (1) SIMULATE n = 100 with low connectance;
            ###  (2) Look only for biggest component (3) and then add
            ### edges as needed 
            
            ###WHEE~ 
             
             
            ### connectance: EXTREMES ARE OK! FOR MODELING. 
             
            
             
             
         }
  }
  
 adj_matrix<-  as_adjacency_matrix(
          deleted_edges_graph,
          type = "both",
          attr = "weight",sparse = F)
  
  
  return(adj_matrix)
}


          
          
  
  
  
  
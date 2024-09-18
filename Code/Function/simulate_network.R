### Simulating the spatial network 
### Subroutines for the main function for simulating the spatial
### network

#' Find the closest values
#'
#' @param vec1_want 
#' @param vec2_have 
#'
#' @return
#' @examples
get_closest_values_vecs <- function(vec1_want, vec2_have){
        index_closest_value <- NULL
        for (value in 1:length(vec1_want)){
                index_closest_value [[value]] <- 
                        which.min(abs(vec1_want[[value]] - vec2_have))
                
        }
        return(do.call(rbind, index_closest_value ))
}

#' Calculate the connectance of a network given the node and edges
#'
#'This calculates the connectance when given a node or edge. To 
#'use this with an igraph, use ecount (number of
#'edges ) and vcount (number of vertices)
#'
#' @param nodes The number of nodes in the network (numeric)
#' @param edges The number of edges in the network (numeric)
#' @return a single numeric connectance value
#' @examples calculate_connectance(nodes = 100, edges = 100)
#' @examples calculate_connectance(nodes = vcount(igraph), edges = ecount(igraph))
calculate_connectance <- function(nodes, edges){
        return(edges/(nodes^2))
}

#' Simulating the x and y coordinates for the spatial network
#'
#'Simulates the xy coordinates of the patches as well as
#'giving the distance matrix across all patches
#'
#' @param seed Default seed is 24601
#' @param max_distance The maximum distance that the xy coordinates lie in
#' @return A list element with the data.frame of xy coordinates and 
#' the distance matrix with kernal
#' @examples simulate_xy_coordinates(seed = 24601, 30)
#'
simulate_xy_coordinates <- function(seed = 24601, max_distance_param){
        set.seed(seed)
        xy <- seq(1, max_distance_param, length.out = 1000)  
        # List of all possible coordinates
        x_coord <- sample(xy, 100, replace = TRUE)  # x-coordinate
        y_coord <- sample(xy, 100, replace = TRUE)  # y-coordinate
        xy_coord <- cbind(x_coord, y_coord)  # xy-coordinates combined
        NegExpDist <- as.matrix(exp(-dist(xy_coord))) # distance matrix with kernel
        
        return(list(xy_coord, NegExpDist))
}

#' Retrieve the component with the most nodes from the simulated network
#'
#'Prunes the original distance matrix and retrieve the biggest component
#'
#' @param xy_list The list from simulate_xy_coordinates that includes both the
#' xy coordinates as well as the matrix (elements being the probability of 
#' dispersal based on the negative exponential dispersal )
#'
#' @return An igraph object
#' @examples retrieve_biggest_component(simulate_xy_coordinates(seed = 24601, 30))
#' 
retrieve_biggest_component <- function(xy_list){
          
         #Convert an adjacency matris to a graph
         Adj_graph <- graph_from_adjacency_matrix(
                  xy_list[[2]], mode = "undirected",
                  diag = FALSE, weighted = TRUE)
          
          # Add latitude and longitude
          V(Adj_graph)$Long <- xy_list[[1]][ ,1]  # x-coordinates
          V(Adj_graph)$Lat <- xy_list[[1]][ ,2]  # y-coordinates
          
          #The number of edges we need for a connectance of 3%
          number_of_edges <- (0.03 * (100^2))
         
          # You choose the top number of edges (from above) and delete
          # everything else
          edge_threshold <- sort(E(Adj_graph)$weight, 
                                 decreasing = TRUE)[number_of_edges]
          
          deleted_edges_graph <- delete_edges(Adj_graph, 
                                              which(E(Adj_graph)$weight < 
                                                            edge_threshold))
          
          # Certain nodes will not be connected to each other and 
          # therefore there will be components
          decomposed_components <- decompose(deleted_edges_graph)
          
          # Count the number of nodes for each component and retrieve
          # the index with the biggest number of nodes 
          biggest_component_length <- which.max(
                  lapply(decomposed_components, function(x) vcount(x))
                  )
          
          #Retrieve our network of interest (biggest component)
          network_bigcomp <- decomposed_components[[biggest_component_length]]
          
          #Rename the vertices
          V(network_bigcomp)$name <- 1:vcount(network_bigcomp)
          
          return(network_bigcomp)
 }

#' Recalculate the distance matrix of the new network 
#'
#' @param network_newcomp The igraph object that is the biggest component from
#' retrieve_biggest_component
#' @return a data.frame that has all the possible edges and their weights in order
#' @examples
recalculate_distance_matrix <- function(network_newcomp) {
        
        # Get the x-y coordinates of the newest component 
        xy_coord_interest <- cbind(
                V(network_newcomp)$Long,
                V(network_newcomp)$Lat)
        
        # Calculate new distance matrices and apply the distance kernal
        DispMat_interest <- as.matrix(exp(-dist(xy_coord_interest)))
        
        #We want to ensure that we're not getting multiedges
        DispMat_interest[lower.tri(DispMat_interest)] <- NA
        
        #The data.frame of already existing edges between the patches
        edgelist_present <- as_edgelist(network_newcomp, names = T)

        ###The new distance matrix
        melted_edge_list <- melt(DispMat_interest, varnames = c("patch1", "patch2"), 
                             value.name = "weight")

        ###Remove self_loop
        melted_edge_list <- na.omit(melted_edge_list[!(melted_edge_list$patch1 == 
                                                    melted_edge_list$patch2),])
        
        
        edges_to_add <- paste0(melted_edge_list$patch1, "-",melted_edge_list$patch2)
        edges_present <- paste0(edgelist_present [,1],"-", edgelist_present [,2])
        
        new_distance <- subset(melted_edge_list,!(edges_to_add %in% edges_present))
        
        new_distance_df <- new_distance[order(new_distance$weight, 
                                              decreasing = TRUE), ]
        
        return(new_distance_df)
}


#' Main function for simulating the spatial network
#'
#'This is the main function that is for simulating the spatial network
#'
#' @param seed Default seed is 24601
#' @param max_distance The maximum distance that the xy coordinates lie in
#' @param specific_connectance A singular connectance value that you are interested
#' in
#' @return A list of igraph 
#' @examples
#' 
simulate_spatial_network <- function(seed, max_distance, specific_connectance =
                                             NA) {
        
        list_xy_coord <- simulate_xy_coordinates(seed, max_distance)
        network_template <- retrieve_biggest_component(list_xy_coord)
        possible_edges_df <- recalculate_distance_matrix(network_template)
        
        adj_list <- NULL
        adj_info_list <- NULL
        
        ### Manually add the first network in
        
        adj_list[[1]] <- network_template
        
        adj_info_list[[1]] <- c(
                
                num_nodes = vcount(network_template),
                num_edges = ecount(network_template),
                
                connectance = calculate_connectance(
                        vcount(network_template),
                        ecount(network_template)
                )
        )
        
        ### For loop time
        for (new_edge in seq(1, nrow(possible_edges_df))) {
                network_template <- add_edges(network_template,
                                 c(possible_edges_df[new_edge, "patch1"], 
                              possible_edges_df[new_edge, "patch2"]),
                       weight = possible_edges_df[new_edge, "weight"]
                )
                
                
                adj_list[[new_edge + 1]] <-network_template
                adj_info_list[[new_edge + 1]] <- c(
                        num_nodes = vcount(network_template),
                        num_edges = ecount(network_template),
                        connectance = round(calculate_connectance(
                                vcount(network_template),
                                ecount(network_template)),5)
                )
        }
        
        adj_info_DF <- as.data.frame(do.call(rbind, adj_info_list))
        
        if (is.na(specific_connectance) == TRUE){
        connectance_interest <- c(adj_info_DF[1,3], seq(0.1, 1,0.1))
        
        interested_network <- get_closest_values_vecs(connectance_interest, 
                                                      adj_info_DF$connectance)
        
        }
        else{
        interested_network <- get_closest_values_vecs(specific_connectance, 
                                                     adj_info_DF$connectance)        
                
        }
        return(adj_list[interested_network])
}


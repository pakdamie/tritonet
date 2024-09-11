###Subroutines for the main function.


#' Calculate the connectance of a network
#'
#' @param nodes The number of nodes 
#' @param edges The number of edges
#'
#' @return a numeric connectance value
#' @export
#'
#' @examples connectance_calculator(nodes = 100,edges = 100)
connectance_calculator <- function(nodes, edges) {
        return(edges / (nodes^2))
}

#' Simulating the x and y coordinates for the spatial network
#'
#' @param seed Default seed is 24601
#' @param max_distance The maximum distance that the xy coordinates lie in
#'
#' @return A list element with the data.frame of xy coordinates and 
#' the distance matrix with kernal
#' @export
#'
#' @examples
#'
simulate_xy_coordinates <- function(seed = 24601, max_distance) {
        set.seed(seed)
        xy <- seq(1, max_distance, length.out = 1000) ### List of all possible coordinates
        x_coord <- sample(xy, 100, replace = TRUE) # x-coordinate
        y_coord <- sample(xy, 100, replace = TRUE) # y-coordinate
        xy_coord <- cbind(x_coord, y_coord) # xy-coordinates combined
        NegExpDist <- as.matrix(exp(-dist(xy_coord))) # distance matrix with kernel
        
        return(list(xy_coord, NegExpDist))
}

tmp <- simulate_xy_coordinates(24601,30)
xy_list <- tmp


#' Retrieve the biggest component from the simulated network
#'
#' @param xy_list The list from simulate_xy_coordinates
#'
#' @return an igraph object
#' @export
#'
#' @examples
retrieve_biggest_component <- function(xy_list) {
        
        Adj_graph <- graph_from_adjacency_matrix(xy_list[[2]],
                                                 mode = "undirected",
                                                 diag = FALSE,
                                                 weighted = TRUE
        )
        
        ## Adding latitude and longitude
        V(Adj_graph)$Long <- xy_list[[1]][,1] # x-coordinates
        V(Adj_graph)$Lat <- xy_list[[1]][,2] # y-coordinates
        
        ###The number of edges we need for a connectance
        number_of_edges <- (0.03 * (100^2))
        
       
        ###You choose the top number of edges (from above) 
        ###and delete everything else 
        deleted_edges_graph <- delete_edges(
                Adj_graph,
                which(E(Adj_graph)$weight < sort(E(Adj_graph)$weight,
                                                 decreasing = T
                )[number_of_edges])
        )
        
        
        decomposed_components <- decompose(deleted_edges_graph)
        
        # Count the number of nodes for each component and then give me
        ### the index for the largest.
        biggest_component_length <- which.max(lapply(
                decomposed_components,
                function(x) {
                        vcount(x)
                }
        ))
        
        ### retrieve our network of interest
        network_bigcomp <- decomposed_components[[biggest_component_length]]
        V(network_bigcomp)$name <- 1:vcount(network_bigcomp)
                    
        return(network_bigcomp)
}



network_newcomp<-        network_bigcomp
#' Recalculate the new network
#'
#' @param network_bigcomp 
#'
#' @return
#' @export
#'
#' @examples
recalculate_distance_matrix <- function(network_newcomp) {
        
        # Get the x-y coordinates of the newest component 
        xy_coord_interest <- cbind(
                V(network_newcomp)$Long,
                V(network_newcomp)$Lat
        )
        
        # Calculate new distance matrices
        DispMat_interest <- as.matrix(exp(-dist(xy_coord_interest)))
        DispMat_interest[lower.tri(DispMat_interest)] <- NA
        #The data.frame of already existing edges between the patches
        edgelist_present <- as_edgelist(network_newcomp, names = T)

        ###The new distance matrix
        melted_edge_list <- melt(DispMat_interest) 
        colnames(melted_edge_list) <- c("patch1", "patch2", "weight")
        
        ###Remove self_loop
        melted_edge_list <- na.omit(melted_edge_list[!(melted_edge_list$patch1 == melted_edge_list$patch2),])
        
        
        edges_to_add <- paste0(melted_edge_list$patch1, "-",melted_edge_list$patch2)
        edges_present <- paste0(edgelist_present [,1],"-", edgelist_present [,2])
        
        new_distance <- subset(
                melted_edge_list,
                !(edges_to_add %in% edges_present)
        )
        
        new_distance_df <- new_distance[order(new_distance$weight, 
                                              decreasing = TRUE), ]
        
        return(new_distance_df)
}


#' Main function for simulating the spatial network
#'
#' @param seed Default seed is 24601
#' @param max_distance The maximum distance that the xy coordinates lie in
#'
#' @return
#' @export
#'
#' @examples
#' 
simulate_spatial_network <- function(seed, max_distance) {
        
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
                
                connectance = connectance_calculator(
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
                        connectance = round(connectance_calculator(
                                vcount(network_template),
                                ecount(network_template)),4)
                )
        }
        
        adj_info_DF <- as.data.frame(do.call(rbind, adj_info_list))
        

        
        
        connectance_interest <- c(adj_info_DF[1,3],seq(0.1, 1,0.1))
        
        
        interested_network <- get_closest_values_vecs(connectance_interest, adj_info_DF $connectance)
        
        
        return( adj_list[interested_network])
}


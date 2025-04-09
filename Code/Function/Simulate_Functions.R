#' Create parameters of interest
#'
#' @returns
#' @export
#'
#' @examples
create_parameters <- function() {
  param_standard <-
    c(
      b_H = 1 / (365 * 70), ## Human mortality rate
      b_P = 1 / 10, # P. Vector birth rate
      b_M = 1 / 10, # S. Vector birth rate
      mu_H = 1 / (365 * 70), ## Human death rate
      mu_V = 1 / 20,
      f_P = 0.25, # Biting rate of the p. vector
      f_M = 0.25, # Biting rate of the s.vector
      theta_P = 0.25, # Transmission probability of p. vector
      theta_M = 0.25 * 0.80, # Transmission probability of s. vector
      theta_H = 0.25, # Transmission probability of human
      gamma = 1 / 90, # Recovery rate of infected human
      c_PM = 3e-4, ## Competition effect of p.vector on s.vector
      c_MP = 3e-6, ## Competition effect of s.vector on p.vector
      c_PP = 4.5e-4, ## Competition effect of p.vector on s.vector
      c_MM = 2.5e-4, ## Competition effect of s.vector on s.vector
      ntime = (365 * 100),
      disturbance_time = (365 * 50),
      delta_T = 1,
      d = 9e-3,
      prop = 0.25,
      mortality_P = 0.25,
      mortality_M = 1
    )




  return(param_standard)
}


#' Create initial states to simulate the model
#'
#' Create the seven matrices for each of the seven compartments. The seven
#' matrices will have the number of rows that depend on time-steps and
#' the number of columns is the number of patches.
#'
#' @param param The parameters that you're choosing
#' @param patch_num The number of patches
#'
#' @return A list of the necessary matrices
#' @export
#'
#' @examples create_initial_states(param_standard, 100)
create_initial_states <- function(
    param,
    patch_num,
    initial_human = 1000,
    initial_vector_s = 2490,
    initial_vector_i = 10,
    infection = "No") {
  compartment_labels <-
    c(
      "HS_mat", "HI_mat", "HR_mat", # Humans (sus/inf/rec)
      "PS_mat", "PI_mat", # Primary (sus/inf)
      "MS_mat", "MI_mat" # Secondary(sus/inf)
    )

  # Using the above compartment label, we then create an empty matrix where
  # the row are the time-steps and the columns are the individual patches
  for (i in 1:length(compartment_labels)) {
    assign(
      compartment_labels[i],
      matrix(0,
        nrow = param["ntime"],
        ncol = patch_num
      )
    )
  }

  initial_vector_i <- ifelse(infection == "No", 0, initial_vector_i)

  # Initial conditions

  # Humans
  HS_mat[1, ] <- rep(initial_human, patch_num)
  HI_mat[1, ] <- rep(0, patch_num)
  HR_mat[1, ] <- rep(0, patch_num)

  # Primary vectors
  PS_mat[1, ] <- rep(initial_vector_s, patch_num)
  PI_mat[1, ] <- rep(initial_vector_i, patch_num)

  # Secondary vectors
  MS_mat[1, ] <- rep(initial_vector_s, patch_num)
  MI_mat[1, ] <- rep(initial_vector_i, patch_num)

  return(list(
    HS_mat, HI_mat, HR_mat,
    PS_mat, PI_mat,
    MS_mat, MI_mat
  ))
}

### Simulating the spatial network

#' Calculate the connectance of a network given the node and edges
#'
#' This calculates the connectance when given a node or edge. To
#' use this with an igraph, use ecount (number of
#' edges ) and vcount (number of vertices)
#'
#' @param nodes The number of nodes in the network (numeric)
#' @param edges The number of edges in the network (numeric)
#' @return a single numeric connectance value
#' @examples calculate_connectance(nodes = 100, edges = 100)
#' @examples calculate_connectance(nodes = vcount(igraph), edges = ecount(igraph))
calculate_connectance <- function(nodes, edges) {
  return(edges / (nodes * (nodes - 1) / 2))
}



#' Simulating the x and y coordinates for the spatial network
#'
#' Simulates the xy coordinates of the patches as well as
#' giving the distance matrix across all patches
#'
#' @param seed Default seed is 24601
#' @param max_distance The maximum distance that the xy coordinates lie in
#' @return A list element with the data.frame of xy coordinates and
#' the distance matrix with kernal
#' @examples simulate_xy_coordinates(seed = 24601, 30)
#'
simulate_xy_coordinates <- function(max_distance, dist_modifer) {
  xy <- seq(1, max_distance, length.out = 1000)

  # List of all possible coordinates
  x_coord <- sample(xy, 100, replace = TRUE) # x-coordinate
  y_coord <- sample(xy, 100, replace = TRUE) # y-coordinate
  xy_coord <- cbind(x_coord, y_coord) # xy-coordinates combined
  NegExpDist <- as.matrix(exp(dist_modifer * -dist(xy_coord))) # distance matrix with kernel

  return(list(xy_coord, NegExpDist))
}

#' Retrieve the component with the most nodes from the simulated network
#'
#' Prunes the original distance matrix and retrieve the biggest component
#'
#' @param xy_list The list from simulate_xy_coordinates that includes both the
#' xy coordinates as well as the matrix (elements being the probability of
#' dispersal based on the negative exponential dispersal )
#'
#' @return An igraph object
#' @examples retrieve_biggest_component(simulate_xy_coordinates(seed = 24601, 30))
#'
retrieve_biggest_component <- function(xy_list) {
  # Convert an adjacency matris to a graph
  Adj_graph <- graph_from_adjacency_matrix(
    xy_list[[2]],
    mode = "undirected",
    diag = FALSE, weighted = TRUE
  )

  # Add latitude and longitude
  V(Adj_graph)$Long <- xy_list[[1]][, 1] # x-coordinates
  V(Adj_graph)$Lat <- xy_list[[1]][, 2] # y-coordinates

  # The number of edges we need for a connectance of 5%
  number_of_edges <- 0.05 * (vcount(Adj_graph) * (vcount(Adj_graph) - 1) / 2)

  # You choose the top number of edges (from above) and delete
  # everything else
  edge_threshold <- sort(E(Adj_graph)$weight, decreasing = TRUE)[number_of_edges]

  deleted_edges_graph <- delete_edges(
    Adj_graph,
    which(E(Adj_graph)$weight <
      edge_threshold)
  )

  # Certain nodes will not be connected to each other and
  # therefore there will be components
  decomposed_components <- decompose(deleted_edges_graph)

  # Count the number of nodes for each component and retrieve
  # the index with the biggest number of nodes
  biggest_component_length <- which.max(
    lapply(decomposed_components, function(x) vcount(x))
  )
  # Retrieve our network of interest (biggest component)
  network_bigcomp <- decomposed_components[[biggest_component_length]]

  # Rename the vertices
  V(network_bigcomp)$name <- 1:vcount(network_bigcomp)

  return(network_bigcomp)
}

#' Recalculate the distance matrix of the new network
#'
#' @param network_newcomp The igraph object that is the biggest component from
#' retrieve_biggest_component
#' @return a data.frame that has all the possible edges and their weights in order
#' @examples
recalculate_distance_matrix <- function(network_newcomp, dist_modifier) {
  # Get the x-y coordinates of the newest component
  xy_coord_interest <- cbind(
    V(network_newcomp)$Long,
    V(network_newcomp)$Lat
  )

  # Calculate new distance matrices and apply the distance kernal
  DispMat_interest <- as.matrix(exp(dist_modifier * -dist(xy_coord_interest)))

  # We want to ensure that we're not getting multiedges
  DispMat_interest[lower.tri(DispMat_interest)] <- NA

  # The data.frame of already existing edges between the patches
  edgelist_present <- as_edgelist(network_newcomp, names = T)

  ### The new distance matrix
  melted_edge_list <-
    melt(DispMat_interest,
      varnames = c("patch1", "patch2"),
      value.name = "weight"
    )

  ### Remove self_loop
  melted_edge_list <- na.omit(melted_edge_list[!(melted_edge_list$patch1 ==
    melted_edge_list$patch2), ])

  edges_to_add <- paste0(melted_edge_list$patch1, "-", melted_edge_list$patch2)
  edges_present <- paste0(edgelist_present[, 1], "-", edgelist_present[, 2])

  new_distance <- subset(melted_edge_list, !(edges_to_add %in% edges_present))

  new_distance_df <- new_distance[order(new_distance$weight,
    decreasing = TRUE
  ), ]

  return(new_distance_df)
}


#' Main function for simulating the spatial network
#'
#' This is the main function that is for simulating the spatial network
#'
#' @param seed Default seed is 24601
#' @param max_distance The maximum distance that the xy coordinates lie in
#' @param specific_connectance A singular connectance value that you are interested
#' in
#' @return A list of igraph
#' @examples
simulate_spatial_network <- function(max_distance, dist_modifier, specific_connectance = NA) {
  # Generate a list of xy coordinates for the network using the seed and max_distance
  list_xy_coord <- simulate_xy_coordinates(max_distance, dist_modifier)

  # Retrieve the largest connected component (network) from the generated coordinates
  network_template <- retrieve_biggest_component(list_xy_coord)

  # Recalculate the distance matrix for potential connections between nodes (patches)
  possible_edges_df <- recalculate_distance_matrix(network_template, dist_modifier)

  # Initialize two lists to store the network structure and related information
  adj_list <- NULL
  adj_info_list <- NULL

  ### Manually add the initial network to the lists
  # Add the initial network (first component) to the adjacency list
  adj_list[[1]] <- network_template

  template_vcount <- vcount(network_template)
  template_ecount <- ecount(network_template)

  # Store network information such as the number of nodes, edges, and connectance
  adj_info_list[[1]] <- c(
    num_nodes = template_vcount, # Number of nodes in the network
    num_edges = template_ecount, # Number of edges in the network

    connectance = calculate_connectance( # Connectance = edges / (nodes * (nodes - 1) / 2)
      nodes = template_vcount,
      edges = template_ecount
    )
  )

  ### Start the loop to add edges to the network for the desired connectance
  # Loop through each potential edge and add it to the network

  # From the template, how many more edges do I need to add in

  for (C in seq(1, length(specific_connectance))) {
    sp_edges_necessary <-
      specific_connectance[C] *
      (template_vcount * (template_vcount - 1) / 2) -
      template_ecount


    # Add an edge between two nodes (patch1 and patch2) with a specific weight
    added_template <- add_edges(network_template,
      c(
        possible_edges_df[1:sp_edges_necessary, "patch1"],
        possible_edges_df[1:sp_edges_necessary, "patch2"]
      ),
      weight = possible_edges_df[1:sp_edges_necessary, "weight"]
    )

    # Update the adjacency list with the new network after adding the edge
    adj_list[[C + 1]] <- added_template

    # Add updated network information to the adj_info_list
    adj_info_list[[C + 1]] <- c(
      num_nodes = vcount(added_template),
      num_edges = ecount(added_template),
      connectance = round(calculate_connectance(
        vcount(added_template),
        ecount(added_template)
      ), 2) # Round connectance to 5 decimal places
    )

    # Convert the adj_info_list to a dataframe for easier analysis
    adj_info_DF <- as.data.frame(do.call(rbind, adj_info_list))
  }


  return(list(adj_list, adj_info_DF))
}



simulate_as_adjacency_matrix <- function(igraph_object) {
  
  adjacency_matrix <- as_adjacency_matrix(
    igraph_object,
    type = "both",
    attr = "weight",
    names = TRUE,
    sparse = FALSE
  )
  summed_prob <- rowSums(adjacency_matrix)
  adjacency_matrix_adj <- sweep(adjacency_matrix, MARGIN = 1, summed_prob, `/`)
  return(adjacency_matrix_adj)
}

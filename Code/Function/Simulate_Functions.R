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
      mu_V = 1/20,
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
      d = 0.05,
      prop = 1,
      mortality_P = 0.25,
      mortality_M = 0.95
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
    initial_vector_s = 1000,
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
  return(edges / (nodes^2))
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
simulate_xy_coordinates <- function(seed = 24601, max_distance_param) {
  set.seed(seed)
  xy <- seq(1, max_distance_param, length.out = 1000)
  # List of all possible coordinates
  x_coord <- sample(xy, 100, replace = TRUE) # x-coordinate
  y_coord <- sample(xy, 100, replace = TRUE) # y-coordinate
  xy_coord <- cbind(x_coord, y_coord) # xy-coordinates combined
  NegExpDist <- as.matrix(exp(5 * -dist(xy_coord))) # distance matrix with kernel

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

  # The number of edges we need for a connectance of 10%
  number_of_edges <- (0.05 * (100^2))

  # You choose the top number of edges (from above) and delete
  # everything else
  edge_threshold <- sort(E(Adj_graph)$weight,
    decreasing = TRUE
  )[number_of_edges]

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
recalculate_distance_matrix <- function(network_newcomp) {
  # Get the x-y coordinates of the newest component
  xy_coord_interest <- cbind(
    V(network_newcomp)$Long,
    V(network_newcomp)$Lat
  )

  # Calculate new distance matrices and apply the distance kernal
  DispMat_interest <- as.matrix(exp(5 * -dist(xy_coord_interest)))

  # We want to ensure that we're not getting multiedges
  DispMat_interest[lower.tri(DispMat_interest)] <- NA

  # The data.frame of already existing edges between the patches
  edgelist_present <- as_edgelist(network_newcomp, names = T)

  ### The new distance matrix
  melted_edge_list <- melt(DispMat_interest,
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
simulate_spatial_network <- function(seed, max_distance, specific_connectance = NA) {
  # Generate a list of xy coordinates for the network using the seed and max_distance
  list_xy_coord <- simulate_xy_coordinates(seed, max_distance)

  # Retrieve the largest connected component (network) from the generated coordinates
  network_template <- retrieve_biggest_component(list_xy_coord)

  # Recalculate the distance matrix for potential connections between nodes (patches)
  possible_edges_df <- recalculate_distance_matrix(network_template)

  # Initialize two lists to store the network structure and related information
  adj_list <- NULL
  adj_info_list <- NULL

  ### Manually add the initial network to the lists

  # Add the initial network (first component) to the adjacency list
  adj_list[[1]] <- network_template

  # Store network information such as the number of nodes, edges, and connectance
  adj_info_list[[1]] <- c(
    num_nodes = vcount(network_template), # Number of nodes in the network
    num_edges = ecount(network_template), # Number of edges in the network

    connectance = calculate_connectance( # Connectance = edges / (nodes * (nodes - 1) / 2)
      nodes = vcount(network_template),
      edges = ecount(network_template)
    )
  )

  ### Start the loop to add edges to the network

  # Loop through each potential edge and add it to the network
  for (new_edge in seq(1, nrow(possible_edges_df))) {
    # Add an edge between two nodes (patch1 and patch2) with a specific weight
    network_template <- add_edges(network_template,
      c(
        possible_edges_df[new_edge, "patch1"],
        possible_edges_df[new_edge, "patch2"]
      ),
      weight = possible_edges_df[new_edge, "weight"]
    )

    # Update the adjacency list with the new network after adding the edge
    adj_list[[new_edge + 1]] <- network_template

    # Add updated network information to the adj_info_list
    adj_info_list[[new_edge + 1]] <- c(
      num_nodes = vcount(network_template),
      num_edges = ecount(network_template),
      connectance = round(calculate_connectance(
        vcount(network_template),
        ecount(network_template)
      ), 5) # Round connectance to 5 decimal places
    )
  }

  # Convert the adj_info_list to a dataframe for easier analysis
  adj_info_DF <- as.data.frame(do.call(rbind, adj_info_list))

  # If a specific connectance is not provided, calculate the range of connectance values
  if (is.na(specific_connectance) == TRUE) {
    connectance_interest <- c(adj_info_DF[1, 3], seq(0.1, 1, 0.1)) # Default connectance range

    # Find the network that is closest to the connectance values of interest
    interested_network <- get_closest_values_vecs(
      connectance_interest,
      adj_info_DF$connectance
    )
  } else {
    # If a specific connectance is provided, find the network closest to that value
    interested_network <- get_closest_values_vecs(
      specific_connectance,
      adj_info_DF$connectance
    )
  }

  # Return the network that matches the closest connectance value
  return(adj_list[interested_network][[1]])
}

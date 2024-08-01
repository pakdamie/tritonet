###These are functions for calculating the adjacency matrix for the modelng


###First, you simulate the adjacency matrix with num_patch being the
###number of patches and connectance being the connectance you would like to 
### see in the map
simulate_adjacency_matrix <- function(seed, num_patch, connectance){
        set.seed(seed)
        num_edges <-  (connectance * (num_patch^2))
        
        ###We create a scale-free network
        new_adjacency <-  erdos.renyi.game(
                num_patch,
                num_edges,
                type = c("gnm"),
                directed = FALSE,
                loops = FALSE
        )
        
        adj_matrix <- as_adjacency_matrix(new_adjacency, sparse = FALSE)
        
}

###Second, you simulate the spatial map with the same amount of num_patch
### what is the max distance of the grid (assume, the villages are laid out
###in an xy cartesian grid)

simulate_spatial_matrix <- function(seed, num_patch, max_distance ){
        set.seed(seed)
        xy <- seq(1,max_distance,length.out=2000)
        x <- sample(xy, num_patch, replace=TRUE)
        y <- sample(xy, num_patch, replace=TRUE)
        lay <- cbind(x,y)
        dst <- as.matrix(exp(-dist(lay)))
        #graph <- graph_from_adjacency_matrix( final_mat, weighted=TRUE, mode='undirected')
        return(dst)
}

###Finally, you combine the adjacency matrix and spatial matrix together
simulate_final_adjacency_matrix <- function(seed, num_patch, connectance, max_distance){
        adj_matrix <- simulate_adjacency_matrix(seed, num_patch,connectance)
        spatial_matrix <- simulate_spatial_matrix(seed, num_patch, max_distance)
        
        return(adj_matrix * spatial_matrix)
}



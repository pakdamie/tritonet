getSparse <- function(g, p=0.25){
        gOut <- igraph::delete_edges(g, 
                                     sample(E(g), round(ecount(g)*p,0), 
                                            prob=1/E(g)$weight)
        )
        return(gOut)
}

# not good at generating spatially-realistic graphs, so use sparsification approaches
# to graphs that make sense spatially instead? 
getSpatial <- function(seed, n, initialOcc=0.5, p=0.10){
        set.seed(seed)
        xy <- seq(1,20,length.out=1000)
        x <- sample(xy, n, replace=TRUE)
        y <- sample(xy, n, replace=TRUE)
        lay <- cbind(x,y)
        dst <- as.matrix(exp(-dist(lay)))
        g <- graph_from_adjacency_matrix(dst, weighted=TRUE, mode='undirected')
        V(g)$Nt <- round(runif(n, 0, 30), 0)
        V(g)$Nt[sample(1:length(V(g)), 
                       round(length(V(g))*(1-initialOcc),0))] <- 0
        g <- getSparse(g, p=p)
        return(g)
}
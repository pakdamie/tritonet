
#' Plot the spatial netowrk graph 
#'
#' @param seed The random seed generator
#' @param num_patch  Number of patches that you want to simulate
#' @param connectance The connectance of the network
#' @param max_distance The maximum distance 
#'
#' @return
#' @export
#'
#' @examples
plot_networkgraph <- function(seed, num_patch, connectance, max_distance){
  
  adjacency_matrix <- simulate_final_adjacency_matrix(seed, num_patch,
                                                      connectance,max_distance)
  
  g9 <- graph_from_adjacency_matrix(adjacency_matrix , weighted=TRUE,
                                    mode="plus", diag=FALSE)
  
  plot(g9)
}

#' Plot the total abundance of the different classes
#'
#' When given the results of discrete_trito_model, this
#' will plot out 7 facets (one for each class) with
#' time on the x-axis and the abundance on the y-axis. Each
#' color reperesents the different patches.
#'
#'
#' @param full_list
#'
#' @return list with the first graph being the complete,full one
#' and the second graph being zoomed in
#' @export
#'
#' @examples
plot_list_groups <- function(full_list) {
  
  identifying_class <- c("HS", "HI", "HR", 
                         "PS", "PI", 
                         "MS", "MI")
  
  data_frame_plot <- lapply(full_list, function(x) {
    data.frame(x, time = seq_len(nrow(x)))
  })
  
  # Assigning the classes so we can combine into one big data.frame
  for (i in seq(1, length(data_frame_plot))) {
    data_frame_plot[[i]]$id <- identifying_class[i]
  }
  
  final_data_frame <- do.call(rbind, data_frame_plot[1:7])
  
  final_melted_frame <- melt(
    final_data_frame,
    id.vars = list("time", "id")
  )
  
  final_melted_frame$id <- factor(
    final_melted_frame$id,
    levels = c(
      "HS", "HI", "HR",
      "PS", "PI", "MS",
      "MI"
    )
  )
  
  eep<- subset(final_melted_frame, final_melted_frame$time > (365 * 50) - 500 &
                 final_melted_frame$time < (365 * 50) + (365*5) )
  
  
  
  full_plot_GG <- ggplot(
    eep,
    aes(x = time, y = (value), color = variable)) +
    geom_line() +
    facet_wrap(~id, scales = "free_y") +
    xlab("Time") +
    ylab("Abundance") +
    theme_classic() +
    scale_color_viridis(discrete = TRUE, option = "turbo") +
    theme(legend.position = "none",
          strip.background = element_blank())
  
  return(full_plot_GG)
}
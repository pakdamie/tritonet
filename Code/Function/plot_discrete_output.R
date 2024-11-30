#' Plot the total abundance of the different classes
#'
#' When given the results of discrete_trito_model, this
#' will plot out 7 facets (one for each class) with
#' time on the x-axis and the abundance on the y-axis. Each
#' color represents the different patches.
#'
#'
#' @param full_list
#'
#' @return a figure
#' @export
#'
#' @examples
plot_list_groups <- function(full_list) {
  
  compartment_labels <- c(
    "HS", "HI", "HR", 
    "PS", "PI", 
    "MS", "MI")
  
  data_frame_plot <- lapply(full_list, function(x) {
    data.frame(x, time = seq_len(nrow(x)))
  })

  # Assigning the classes so we can combine into one big data.frame
  for (i in seq(1, length(data_frame_plot))) {
    data_frame_plot[[i]]$id <- compartment_labels[i]
  }

  final_data_frame <- do.call(rbind, data_frame_plot)

  final_melted_frame <- melt(
    final_data_frame,
    id.vars = list("time", "id")
  )

  final_melted_frame$id <- factor(
    final_melted_frame$id,
    levels = compartment_labels
  )

  full_plot_GG <- ggplot(
    final_melted_frame,
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

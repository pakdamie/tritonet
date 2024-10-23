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
    
  #Convert this to data.frames
  data_frame_plot <- lapply(
    full_list, function(x) data.frame(x))
  
  data_frame_plot <- lapply(data_frame_plot, function(x) {
    x$time <- seq(1, nrow(x)) 
    return(x)}
  )
  
  identifying_class <- c("HS", "HI", "HR", "PS", "PI", "SS", "SI")
  
  #Assigning the classes so we can combine into one big data.frame
  for (i in seq(1,length(data_frame_plot))){
      data_frame_plot[[i]]$id <- identifying_class[i]
     
  }
  
  final_data_frame <- do.call(rbind, data_frame_plot)
  
  final_melted_frame <- melt(
    final_data_frame, 
    id.vars = list('time','id'))
        
  final_melted_frame$id <- factor(
    final_melted_frame$id, levels = c("HS","HI","HR",
                                     "PS","PI", "SS",
                                     "SI"))
  
  
  full_plot_GG <- ggplot(final_melted_frame, 
    aes(x = time, y = (value), color = variable))+
    geom_line() +
    facet_wrap(~id, scales = 'free_y') + 
    xlab("Time")+
    ylab("Abundance")+
    theme_classic()+
    scale_color_viridis(discrete = TRUE, option = 'turbo')+
    theme(legend.position = 'none')  
           
  zoomed_plot_GG <-
      ggplot(final_melted_frame, 
        aes(x = time, y = (value), color = variable))+
      geom_line(alpha= 0.5) +
      xlim(480,1000)+
      facet_wrap(~id) + 
      xlab("Time")+
      ylab("Abundance")+
      theme_classic()+
      scale_color_viridis(discrete = TRUE, option = 'turbo')+
      theme(legend.position = 'none')  
  
  return(list(full_plot_GG,zoomed_plot_GG ))
}



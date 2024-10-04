plot_list_groups <- function(list){
    data_frame_plot <-  data.frame(list)
    data_frame_plot$time <- seq(1, nrow(data_frame_plot))
    
    melted_frame <- melt(data_frame_plot, id.vars = 'time')
        
    ggplot(melted_frame, aes(x = time, y= value, color = variable))+
            geom_line() +
            xlab("Time")+
            ylab("Abundance")+
            theme_classic()+
            scale_color_viridis(discrete = TRUE)+
            theme(legend.position = 'none')  
           
        
}


plot_R0_groups <- function(list){
        
       test_dataframe <- lapply(list, function(x) as.data.frame(x))
     
       for(k in seq(1,length(test_dataframe ))){
               test_dataframe [[k]]$patch <- seq(1,nrow(test_dataframe [[k]]))
               test_dataframe [[k]]$time <- k
               
       }

       full_R0_dataframe <- do.call(rbind,test_dataframe)
       
       ggplot(full_R0_dataframe, aes(x = time, y= (V1), color = as.factor(patch), 
                                     group = as.factor(patch)))+
               geom_point(size =1) +
               scale_color_viridis(discrete = TRUE,option ='inferno') +
               theme_classic() +
               xlab("Time") + 
               ylab("R effective")+
               theme(legend.position = 'none')
}
list <- test_function_L[[1]]


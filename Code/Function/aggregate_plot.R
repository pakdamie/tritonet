aggregator_list <- function(big_list = HIGH_NETWORK_LIST){
        
        
  for (network in seq(1,100)){
          for (entry in 1:9){
                  big_list [[network]][[entry]]$networkid <- network
          }       
  }

        
  full_list <- unlist(big_list, recursive = FALSE)
  
  dominance_region <- lapply(full_list,function(x) cbind(calculate_dominance_region(x),
                             x[,c('frequency','coverage',"networkid")]))
        
  dominance_region  <- lapply(dominance_region ,function(x) {
          x$dominant_species <-ifelse(x$prop >=0.5, "P", "S")
                                               return(x)})
  
  dominance_region_f <- do.call(rbind, dominance_region)
  
  a<-ggplot(dominance_region_f, aes(x = time, y= prop,color = networkid, group =
                                         networkid))+
          scale_color_viridis()+
          geom_point(alpha = 0.5)+ 
          geom_line()+
          geom_hline(yintercept = 0.50)+
          theme_bw()+
          facet_grid(frequency~coverage)
  
  
  time_dom_region <-  do.call(rbind,lapply(dominance_region, function(x){
                                 table(x$dominant_species)["P"]/sum(table(x))}))
  
  freq_coverage <- unique(do.call(rbind, lapply(dominance_region, function(x)
                           x[,c("frequency","coverage","networkid")])))
  
  aggregated_df <- cbind(time_dom_region, freq_coverage)
  
  aggregated_df$P[is.na(aggregated_df$P)== TRUE] <- 0
   
  
  
  
  gg_plot <- ggplot(aggregated_df, aes(x = as.factor(frequency), y= P))+
          geom_boxplot(aes(fill= as.factor(coverage)),position="dodge", width = 0.4) + 
          scale_fill_viridis(option = 'mako',
                             discrete = TRUE,
                             name = "Coverage")+theme_classic()+
          theme(axis.text = element_text(size = 12), 
                axis.title = element_text(size = 14))+
          ggtitle("Connectance - 50%")+
          xlab("Frequency")+
          ylab("Proportion of time dominant");
  a+gg_plot
  return(gg_plot)
}  
  
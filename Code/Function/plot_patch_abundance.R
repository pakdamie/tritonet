
#' Reformat the list of desolve outputs for plotting (vectors)
#'
#' @param list The desolve outputs that are in a list
#'
#' @return A data.frame with the different abundance of PS,PI,SS, and SI with time
#' as well as the patch numbers.
#' @export
#'
#' @examples
reformat_abundance_patch_vector <- function(list){
        
      abundance_matrices <-  calculate_vector_abundance_patch(list) 
        
      melted_abundance_matrices<- lapply(abundance_matrices, 
                                         function(x) melt(as.data.frame(x), id.vars = "time"))
      
      melted_abundance_matrices <- lapply(melted_abundance_matrices, 
                                          function(x) 
                                        cbind.data.frame(x, 
                                                         idvariable = gsub("[^a-zA-Z]", "", x$variable),
                                                         patch_num = parse_number(as.character(x$variable))))
     
      full_abundance_vectors <- do.call(rbind,melted_abundance_matrices )
      
      
      ###This is for the creation of data.frame that has the total abundance
      ### of all primary (regardless of status) and secondary
      
      PSI_only_df <- subset(full_abundance_vectors ,  
                            full_abundance_vectors $idvariable %in% c("PS", "PI"))
      
      agg_PSI_only_df <- aggregate(PSI_only_df$value, list(time =  PSI_only_df$time,
                                                            patch_num =PSI_only_df$patch_num),
                                    'sum')
      
      agg_PSI_only_df$id <- "Primary"
      
      SSI_only_df <- subset(full_abundance_vectors ,  
                            full_abundance_vectors $idvariable %in% c("SS", "SI"))
      
       
      agg_SSI_only_df  <- aggregate( SSI_only_df$value, list(time =  SSI_only_df$time,
                                                             patch_num = SSI_only_df$patch_num),
                                     'sum')
      agg_SSI_only_df$id <- "Secondary"
      
      agg_full_df <- rbind(agg_PSI_only_df,   agg_SSI_only_df )
      
      

       return(list(full_abundance_vectors, agg_full_df ))
}

#' Reformat the list of desolve outputs for plotting (hosts)
#'
#' @param list 
#'
#' @return
#' @export
#'
#' @examples
reformat_abundance_patch_host <- function(list){
        
        abundance_matrices <-  calculate_host_abundance_patch(list) 
        
        melted_abundance_matrices<- lapply(abundance_matrices, 
                                           function(x) melt(as.data.frame(x), id.vars = "time"))
        
        melted_abundance_matrices <- lapply(melted_abundance_matrices, 
                                            function(x) 
                                                    cbind.data.frame(x, 
                                                                     idvariable = gsub("[^a-zA-Z]", "", x$variable),
                                                                     patch_num = parse_number(as.character(x$variable))))
        
        full_abundance_host <- do.call(rbind,melted_abundance_matrices)
        
    
        
        return(full_abundance_host)
}


#' Plot vector abundance 
#'
#' @param list 
#' @param aggregate 
#'
#' @return
#' @export
#'
#' @examples
plot_abundance_patch_time <- function(list, aggregate = "no"){
      
      gg_df <-  reformat_abundance_patch_vector(list)[[1]]
      agg_full_df <- reformat_abundance_patch_vector(list)[[2]]
              
      gg_abundance_full <- ggplot(agg_full_df,    
                                  aes(x = time, y= log10(x+1), 
                                      color = patch_num,
                                      group =  patch_num))+
              geom_line(alpha = 0.6)+
              facet_wrap(~ id,ncol = 1)+
              xlab("Time")+
              ylab("Abundance ")+
              scale_color_viridis(option = 'viridis') + 
              theme_minimal()
              
              
      gg_abundance_facet <- ggplot(gg_df, 
                             aes(x = time, y= log10(value+1), 
                                 color = patch_num,
                                 group =  variable))+
              geom_line()+
              facet_wrap(~ idvariable)+
              xlab("Time")+
              ylab("Abundance ")+
              scale_color_viridis(option = 'viridis') + 
              theme_minimal()
         
        
      if(aggregate == "no"){
      
       return( gg_abundance_full + gg_abundance_facet + plot_layout(guides = 'collect'))
      }
      else{
              return(gg_abundance_full)
      }
}

        
plot_abundance_patch_time_host <- function(list){
        
        gg_df <-  reformat_abundance_patch_host(list)
        
        gg_abundance_full <- ggplot(subset(gg_df,gg_df$idvariable!= "HS"),    
                                    aes(x = time, y = log10(value+1), 
                                        color = patch_num,
                                        group =  patch_num))+
                geom_line(alpha = 0.6)+
                facet_wrap(~ idvariable,ncol = 3)+
                xlab("Time")+
                ylab("Abundance ")+
                scale_color_viridis(option = 'viridis') + 
                theme_minimal()
        
       return( gg_abundance_full)
}

calculate_max_CV_RE <-  function(list){

  test_dataframe <- lapply(list, function(x) as.data.frame(x))
        
        
  for(k in seq(1,length(test_dataframe ))){
                test_dataframe [[k]]$patch <- seq(1,nrow(test_dataframe [[k]]))
                test_dataframe [[k]]$time <- k
                
        }
   full_list <- do.call(rbind, lapply(test_dataframe , function (x) {
     data.frame(max_RE = max(x$V1), 
     CV_RE = sd(x$V1)/mean(x$V1))}))
        
   
   full_list_post300 <- full_list[300:1000,]
     
   return(data.frame(max = max(
           full_list_post300 $max_RE,na.rm = TRUE),
         cv = max(full_list_post300 $CV_RE, na.rm = TRUE)))
        
        
        
        
        
}

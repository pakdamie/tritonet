#' Calculate the maximum R-effective number
#'
#' @param list (The output)
#'
#' @return
#' @export
#'
#' @examples
calculate_max_CV_RE <-  function(list){

  #Turn into data.frame so we can add patch number and time
  test_dataframe <- lapply(list, function(x) as.data.frame(x))
        
        
  for(k in 1:length(test_dataframe)) {
    test_dataframe [[k]]$patch <- seq(1,nrow(test_dataframe [[k]]))
    test_dataframe [[k]]$time <- k
  }
  
   full_list <- do.call(
     rbind, 
     lapply(
       test_dataframe , function (x) {
         data.frame(
           max_RE = max(x$V1), CV_RE = sd(x$V1)/mean(x$V1))}))
        
   #R-effective/R0 will always be the greatest at the beginning
   #so just get it here
   full_list_post300 <- full_list[300:1000,]
     
   return(data.frame(max = max(
           full_list_post300 $max_RE,na.rm = TRUE),
         cv = max(full_list_post300 $CV_RE, na.rm = TRUE)))
        
}

calculate_CV_Reff_patch <- function(list_reffective){
        
        test_dataframe <- lapply(list_reffective, 
                                 function(x) as.data.frame(x))
        
        for(k in seq(1,length(test_dataframe ))){
                test_dataframe [[k]]$patch <- seq(1,nrow(test_dataframe [[k]]))
                test_dataframe [[k]]$time <- k
                
        }
        
        full_R0_dataframe <- do.call(rbind,test_dataframe)
        
        split_patch <- split(full_R0_dataframe, 
                             list(full_R0_dataframe$time))
        
        CV_R_Effective <- lapply(split_patch, function(x) sd(x$V1)/mean(x$V1))

        full_CV <- data.frame(time = 1:999,
                   CV = do.call(rbind, CV_R_Effective))
        

        
        return(full_CV)
        
}


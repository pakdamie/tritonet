c_PP_standard <- 4.5e-4 ## Competition effect of p.vector on s.vector
c_MM_standard <- 2.5e-4 ## Competition effect of s.vector on p.vector
c_PM_standard <- 2.5e-4 ## Comeptition effect of primary on s vector
c_MP_standard <- 4.5e-4  ## Comeptition effect of secondary on p.vector

C_abundance<- function(c_PP_standard , c_MM_standard ,c_PM_standard,c_MP_standard ){
        
 param_standard <- get_parameters("standard")       
        
 param_standard["c_PP"] <- c_PP_standard 
 param_standard["c_MM"] <- c_MM_standard 
 param_standard["c_PM"] <- c_PM_standard 
 param_standard["c_MP"] <- c_MP_standard         
 param_standard["mortality_P"] <- c_MP_standard         
 
output <- Simulate_Model_Output(
   parameter =  param_standard,
   infection_start = "No",
   variable_interest = NA,
   vector_value = NA) |>
 Calculate_Human_Reff_Expanded(param =  param_standard)
 
return(output)     
        
}


tmp <- C_abundance(4.5e-4,4.5e-4,4.5e-4 * 0.1, 4.5e-4 * 0.1)
tmp_sub<- subset(tmp, tmp$time > 9000 & tmp$time < 9125 + 1200)

ggplot(tmp_sub, aes(x = time, y =NP + NM)) + geom_path()


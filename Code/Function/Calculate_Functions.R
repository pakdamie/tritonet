#' Title
#'
#' @param parameters 
#' @param lists 
#' @param disturbance_time 
#'
#' @returns
#' @export
#'
#' @examples
calculate_R_effective_discrete_net <- function(parameters, lists, disturbance_time) {
  
  # Extract parameters
  theta_P <- parameters["theta_P"]
  theta_M <- parameters["theta_M"]
  theta_H <- parameters["theta_H"]
  gamma <- parameters["gamma"]
  mu_H <- parameters["mu_H"]
  mu_V <- parameters["mu_V"]
  f_P <- parameters["f_P"]
  f_M <- parameters["f_M"]
  c_MP <- parameters["c_MP"]
  c_PM <- parameters["c_PM"]
  c_MM <- parameters["c_MM"]
  c_PP <- parameters["c_PP"]
  
  # Extract lists
  HS <- lists[[1]]
 
  PS <- lists[[4]]
  PI <- lists[[5]]
  
  MS <- lists[[6]]
  MI <- lists[[7]]
  
  patch_num <- ncol(HS)
  
  NH <- 1000
  NP <- PS + PI
  NM <- MS + MI
  
  HtoP <- (f_P * theta_H)/(gamma + mu_H) * (PS/NH)
  PtoH <- (f_P * theta_P)/(mu_V) * (HS/NH)

  HtoM <- (f_M * theta_H)/(gamma + mu_H) * (MS/NH)
  MtoH <- (f_M * theta_M)/(mu_V) * (HS/NH)
  
  Patch_RE <- (HtoP * PtoH) + (HtoM * MtoH)
  
  return(cbind.data.frame(time = 1:nrow(HS), Patch_RE, total_RE = rowSums(Patch_RE)))

}


check_stability <- function(full_list){
  
  
  
  
  
  
}


calculate_CV_RE <- function(df){
  
  
  apply()
  
  
  
  
  
  
}









namer_chosen_compartments <- function(chosen_patch) {
  
  compartments <- c("PS", "PI", #primary susceptible and infected
                    "SS", "SI") #secondary susceptible and infected
  
  # Generate the names for each compartment
  result <- lapply(compartments, function(prefix) {
    apply(chosen_patch, 2, function(x) paste0(prefix, x), 
          simplify = TRUE)
  })
  
  
  return(do.call(rbind,result))
}










#' Find the closest values
#'
#' @param vec1_want 
#' @param vec2_have 
#'
#' @return
#' @examples
get_closest_values_vecs <- function(vec1_want, vec2_have){
  index_closest_value <- NULL
  
  for (value in 1:length(vec1_want)){
    index_closest_value [[value]] <- 
      which.min(abs(vec1_want[[value]] - vec2_have))
  }
  return(do.call(rbind, index_closest_value ))
}

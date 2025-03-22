#' Calculate the reproductive number (human) derived analytically 
#' assuming that the susceptible human is constant and only the 
#' susceptible vector population is changing. This is the main
#' way that the RE was calculated in the main manuscript
#'
#' @param List_x (list) The raw model output that should include 7 elements
#' @param param (data.frame) The data.frame with the parameter
#'
#' @returns (data.frame) A large data.frame with the RE as well as the 
#' different "parts" that describe the RE
#' @export
#'
#' @examples Calculate_Human_Reff_Expanded(Model_List, param_standard)

Calculate_Human_Reff_Expanded <- function(List_x, param) {
  if (length(List_x) != 7 | class(List_x) != "list") {
    stop("This is not the correct format:
         x`are you sure this is the model output?")
  }

  # Parameters
  theta_H <- param[["theta_H"]] # Transmission probability of human
  theta_P <- param[["theta_P"]] # Transmission probability of primary vector
  theta_M <- param[["theta_M"]] # Transmission probability of secondary vector
  f_P <- param[["f_P"]] # Primary biting rate
  f_M <- param[["f_M"]] # Secondary biting rate
  gamma <- param[["gamma"]] # Recovery rate
  mu_H <- param[["mu_H"]] # Human mortality rate
  mu_V <- param[["mu_V"]] # Vector mortality rate
  c_MP <- param[["c_MP"]] # Competition of secondary on primary
  c_PM <- param[["c_PM"]] # Competition of primary on secondary
  c_MM <- param[["c_MM"]] # Competition of secondary on secondary
  c_PP <- param[["c_PP"]] # Competition of primary on primary
  ntime <- param[["ntime"]] # How long does the simulation run for?

  # We keep the number of susceptible human as being consistent 
  HS <- 1000 #number of susceptible human
  NH <- 1000 #total human 

  PS <- List_x[[4]] # primary susceptible
  PI <- List_x[[5]] # primary infected

  MS <- List_x[[6]] # secondary susceptible
  MI <- List_x[[7]] # secondary infected

  NP <- PS + PI # Total primary vectors
  NM <- MS + MI # Total secondary vectors

  # Related to the average dwell time of the human
  human_wait_time <- gamma + mu_H

  #Primary vector contribution to the RE
  Primary_HP <- (theta_H * f_P * (PS / NH)) * (1 / human_wait_time)
  Primary_PH <- (theta_P * f_P * (HS / NH)) * (1 / mu_V)

  #Secondary contribution to the RE
  Secondary_HM <- (theta_H * f_M * (MS / NH)) * (1 / human_wait_time)
  Secondary_MH <- (theta_M * f_M * (HS / NH)) * (1 / mu_V)

  #The total RE
  RE <- (Primary_HP * Primary_PH) + (Secondary_HM * Secondary_MH)


  RE_DF <- cbind.data.frame(
    time = seq(1, ntime), #Time
    RE = RE, #RE
    NP = NP, #Total primary vector 
    NM = NM, #Total secondary vector 
    PtoH = (Primary_HP * Primary_PH), # P.vector contribution
    MtoH = (Secondary_HM * Secondary_MH), #S.vector contribution
    Primary_HP = Primary_HP, #Human to primary vector
    Primary_PH = Primary_PH, #P.vector to human 
    Secondary_HM = Secondary_HM, #Human to secondary vector
    Secondary_MH = Secondary_MH #S.vector to human
  )
  return(RE_DF)
}

Calculate_Human_Reff_Static<- function(P,M, param) {
 
  
  # Parameters
  theta_H <- param[["theta_H"]] # Transmission probability of human
  theta_P <- param[["theta_P"]] # Transmission probability of primary vector
  theta_M <- param[["theta_M"]] # Transmission probability of secondary vector
  f_P <- param[["f_P"]] # Primary biting rate
  f_M <- param[["f_M"]] # Secondary biting rate
  gamma <- param[["gamma"]] # Recovery rate
  mu_H <- param[["mu_H"]] # Human mortality rate
  mu_V <- param[["mu_V"]] # Vector mortality rate
  c_MP <- param[["c_MP"]] # Competition of secondary on primary
  c_PM <- param[["c_PM"]] # Competition of primary on secondary
  c_MM <- param[["c_MM"]] # Competition of secondary on secondary
  c_PP <- param[["c_PP"]] # Competition of primary on primary
  ntime <- param[["ntime"]] # How long does the simulation run for?
  
  # We keep the number of susceptible human as being consistent 
  HS <- 1000 #number of susceptible human
  NH <- 1000 #total human 
  
  NP = PS = P
  
  NM =MS = M
  # Related to the average dwell time of the human
  human_wait_time <- gamma + mu_H
  
  #Primary vector contribution to the RE
  Primary_HP <- (theta_H * f_P * (PS / NH)) * (1 / human_wait_time)
  Primary_PH <- (theta_P * f_P * (HS / NH)) * (1 / mu_V)
  
  #Secondary contribution to the RE
  Secondary_HM <- (theta_H * f_M * (MS / NH)) * (1 / human_wait_time)
  Secondary_MH <- (theta_M * f_M * (HS / NH)) * (1 / mu_V)
  
  #The total RE
  RE <- (Primary_HP * Primary_PH) + (Secondary_HM * Secondary_MH)
  
  
  RE_DF <- cbind.data.frame(
    RE = RE, #RE
    NP = NP, #Total primary vector 
    NM = NM #Total secondary vector 

  )
  return(RE_DF)
}






#' Calculate the total vector abundance for the maximum RE AND the
#' RE for the maximum vector abundance.
#'
#' Given a data.frame, we calculate by group (decided by the splitting_factor)
#' what is the maximum RE and the corresponding vector abundance as well as
#' what is the maximum vector abundance and the corresponding RE? Function is
#' to tackle the question of is the max vector abundance always the max RE?
#'
#' @param df (data.frame) The data.frame where the RE is calculated
#' @param splitting_factor What variable are we comparing groups by?
#'
#' @returns (list) where the first element is the maximum RE, and the second
#' element is the maximum vector abundance and the corresponding RE
#' @export
#'
#' @examples
Calculate_max_RE_DF <- function(df, splitting_factor) {
  split_df <- split(df, df[, splitting_factor])

  # Maximum RE and the corresponding vector abundance
  max_RE_DF <- do.call(
    rbind,
    lapply(split_df, function(x) x[which.max(x$RE), ])
  )

  min_RE_DF <- do.call(
    rbind,
    lapply(split_df, function(x) x[which.min(x$RE), ])
  )

  # Maximum vector abundance and corresponding RE
  max_Vab_DF <- do.call(
    rbind,
    lapply(split_df, function(x) x[which.max(x$NP + x$NM), ])
  )

  return(list(max_RE_DF, min_RE_DF, max_Vab_DF))
}

#' Calculate the increase above the baseline (equilibrium)
#' )
#' Given a data.frame, we calculate by group (decided by the splitting_factor)
#' what is the maximum RE and the corresponding vector abundance as well as
#' what is the maximum vector abundance and the corresponding RE? Function is
#' to tackle the question of is the max vector abundance always the max RE?
#'
#' @param model_output (list) Raw model output with abundance
#' @param parameter_list (list) Collection of parameters (data.frame)
#' @param parameter_list (data.frame) The expanded.grid form
#' @param temp (character)
#'
#' @returns (data.frame) 
#'
#' @examples
Calculate_change_baseline <- function(model_output, 
                                      parameter_list, 
                                      expand_param, 
                                      temp) {
  
  param_dsturb <- parameter_list[[1]][["disturbance_time"]]
  
  
  
  RE_list <- NULL
  for (i in seq(parameter_list)) {
    RE_tmp <- Calculate_Human_Reff_Expanded(model_output[[i]], 
                                            parameter_list[[i]])

    eq_RE <-   RE_tmp[param_dsturb - 1, ]$RE
    eq_NM <-   RE_tmp[param_dsturb -1, ]$NM
    eq_NP <-   RE_tmp[param_dsturb -1, ]$NP
    eq_mtoH <- RE_tmp[param_dsturb -1, ]$MtoH
    
    
  
    
    post_df <- na.omit(subset(RE_tmp,
                              RE_tmp$time > param_dsturb - 100 
                              & RE_tmp$time < param_dsturb + 500))

    
    RE_max <- cbind(expand_param[i, ],
      max_NP = max(post_df$NP) - eq_NP,
      max_NM = max(post_df$NM) - eq_NM,
      max_NV = max(post_df$NM + post_df$NP) - (eq_NM + eq_NP),
      RE = max(post_df$RE) - eq_RE,
      max_RE = max(post_df$RE),
      max_mtoH = max(post_df$MtoH/post_df$RE)
    )
    
    RE_max_temp <- cbind(expand_param[i, ],
      time = post_df$time,
      max_NP = post_df$NP ,
      max_NM = post_df$NM,
      max_NV =post_df$NM + post_df$NP,
      RE = post_df$RE,
      max_mtoH = post_df$MtoH/post_df$RE
    )
    
    RE_list[[i]] <- switch(temp,
                           "No" =  RE_max ,
                           "Yes" = RE_max_temp)
  }
  return(do.call(rbind, RE_list))
}

Isocline_Primary <- function(c_PM_input) {
  ab <- seq(0, 1200)
  b_V <- param_standard["b_P"]
  mu_V <- param_standard["mu_V"]
  
  c_PP <- param_standard["c_PP"]
  c_MP <- param_standard["c_MP"]
  c_MM <- param_standard["c_MM"]
  c_PM <- c_PM_input
  
  isocline_P_P0 <- ((1 - mu_V / b_V)) / c_MP
  isocline_P_M0 <- ((1 - mu_V / b_V)) / c_PP
  isocline_M_P0 <- ((1 - mu_V / b_V)) / c_MM
  isocline_M_M0 <- ((1 - mu_V / b_V)) / c_PM
  
  # Corrected isocline placements
  isocline_P <- data.frame(x = 0, y =isocline_P_P0 , xend =   isocline_P_M0, yend = 0, id = "P")
  isocline_M <- data.frame(x = 0, y = isocline_M_P0, xend = isocline_M_M0 , yend = 0, id = "M")
  
  return(list(isocline_P, isocline_M))
}


Phaseplot_Primary <- function(c_PM_input) {
  P <- seq(0,1100,100)
  M <- seq(0,1800,100)
  expanded_PM <- expand.grid(P,M)
  b_V <- param_standard["b_P"]
  mu_V <- param_standard["mu_V"]
  
  
  
  c_PP <- param_standard["c_PP"]
  c_MP <- param_standard["c_MP"]
  c_MM <- param_standard["c_MM"]
  c_PM <- c_PM_input
  
  
  slope_P <- b_V *  expanded_PM $Var1 * (1- c_PP * expanded_PM $Var1 - c_MP * expanded_PM $Var2) - mu_V *expanded_PM $Var1
  slope_M <- b_V *  expanded_PM $Var2  * (1- c_MM *  expanded_PM $Var2  - c_PM *  expanded_PM $Var1 ) - mu_V * expanded_PM $Var2 
  
  
  RE<- Calculate_Human_Reff_Static(expanded_PM $Var1,expanded_PM $Var2, get_parameters("standard"))
  
  
  
  
  abundance_slope<- cbind(expanded_PM, slope_P, slope_M, RE = RE$RE)
  
  return(abundance_slope)
}


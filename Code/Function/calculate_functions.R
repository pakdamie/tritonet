#' Calculate R_effective (human) analytically
#'
#' Calculate the effective reproductive number with the assumption that any 
#' new infection is ONLY from the human (check the supplementary material).
#' 
#' @param List_x A list that come from the model output. Should have length
#' of 7
#' @param param A data.frame of parameters to be used for calculating the model
#'
#' @return A data.frame that includes the RE, the total number of primary (NP)
#' and secondary vectors (NM)
#' 
#' @export
#'
#' @examples Calculate_Human_REff(model_output, param_standard)
#' 
Calculate_Human_REff <- function(List_x, param) {
        
  # Parameters        
  theta_H <- param["theta_H"] #Transmission probability of human
  theta_P <- param["theta_P"] #Transmission probability of primary vector
  theta_M <- param["theta_M"] #Transmission probability of secondary vector
  f_P <- param["f_P"] # Primary biting rate
  f_M <- param["f_M"] # Secondary biting rate
  gamma <- param["gamma"] # Recovery rate
  mu_H <- param["mu_H"] # Human mortality rate
  c_MP <- param["c_MP"] #Competition of secondary on primary
  c_PM <- param["c_PM"] #Competition of primary on secondary
  c_MM <- param["c_MM"] #Competition of secondary on secondary
  c_PP <- param["c_PP"] #Competition of primary on primary
  ntime <- param["ntime"] #How long does the simulation run for?
  
  # States
  HS <- List_x[[1]]
  HI <- List_x[[2]]
  HR <- List_x[[3]]
  
  PS <- List_x[[4]]
  PI <- List_x[[5]]
  
  MS <- List_x[[6]]
  MI <- List_x[[7]]
  
  NH <- HS + HI + HR # Total humans
  NP <- PS + PI # Total primary vectors
  NM <- MS + MI # Total secondary vectors
  
  Psi_P <- (c_MP * NM) + c_PP * (PS + (2 * PI))
  Psi_M <- (c_PM * NP) + c_MM * (MS + (2 * MI))
  Psi_C <- (c_MP * PI) * (c_PM * MI)
  Gamma <- (Psi_P * Psi_M) - Psi_C
  
  ### Waiting time denominator 
  wait_time <- Gamma * (gamma + mu_H) * NH^2
  
  Primary <- HS * (theta_H * f_P * PS) * 
             ((Psi_M * (theta_P * f_P) - (c_PM * MI * (theta_M * f_M))))
  Secondary <- HS*  (theta_H * f_M * MS) * 
             ((Psi_P * (theta_M * f_M) - (c_MP * PI * (theta_P * f_P))))
  
  RE <- (Primary / wait_time) + (Secondary / wait_time) 
  
  RE_DF <- cbind.data.frame(
    time = seq(1, ntime),
    RE = RE,
    NP = NP,
    NM = NM,
    PS = PS,
    MS = MS,
    HS = HS,
    PtoH = Primary/wait_time,
    MtoH = Secondary/wait_time,
    wait_time = wait_time)
  
  return(RE_DF)
}

#' Calculate R_0 (human) analytically
#'
#' Function is different than `Calculate_Human_Reff`
#'
#' @param df_expand  
#' @param param 
#'
#' @returns
#' @export
#'
#' @examples
Calculate_Human_Reff_Expanded <- function(df_expand, param){
        
  # Parameters        
  theta_H <- param["theta_H"] #Transmission probability of human
  theta_P <- param["theta_P"] #Transmission probability of primary vector
  theta_M <- param["theta_M"] #Transmission probability of secondary vector
  f_P <- param["f_P"] # Primary biting rate
  f_M <- param["f_M"] # Secondary biting rate
  gamma <- param["gamma"] # Recovery rate
  mu_H <- param["mu_H"] # Human mortality rate
  c_MP <- param["c_MP"] #Competition of secondary on primary
  c_PM <- param["c_PM"] #Competition of primary on secondary
  c_MM <- param["c_MM"] #Competition of secondary on secondary
  c_PP <- param["c_PP"] #Competition of primary on primary

  HS <- 1000
  NH <- 1000
  
  NP <- df_expand["NP"]
  PI <- NP * 0.10  
  PS <- NP - PI
  
  NM <- df_expand["NM"]
  MI <- NM * 0.10  
  MS <- NM - MI
  
  
  Psi_P <- (c_MP * NM) + c_PP * (PS + (2 * PI))
  Psi_M <- (c_PM * NP) + c_MM * (MS + (2 * MI))
  Psi_C <- (c_MP * PI) * (c_PM * MI)
  Gamma <- (Psi_P * Psi_M) - Psi_C
  
  wait_time <- Gamma * (gamma + mu_H) * NH^2
  
  Primary <-   HS * (theta_H * f_P * PS) * ((Psi_M * (theta_P * f_P) - 
               (c_PM * MI * (theta_M * f_M))))
  Secondary <- HS*  (theta_H * f_M * MS) * ((Psi_P * (theta_M * f_M) - 
               (c_MP * PI * (theta_P * f_P))))
  
  
  RE <- (Primary/wait_time) + (Secondary/wait_time) 
  
  RE_DF <- cbind.data.frame(
          RE = RE,
          NP = NP,
          NM = NM,
          PS = PS,
          MS = MS,
          HS = HS,
          PtoH = Primary/wait_time,
          MtoH = Secondary/wait_time
  )
  
  colnames(RE_DF) <- c("RE", "NP","NM", "PS","MS","HS","PtoH","MtoH")
  
  return(RE_DF)
}


#' Calculate the vector abundance for the maximum RE or the
#' RE for the maximum vector abundance
#'
#' @param df 
#' @param splitting_factor 
#'
#' @returns A list where the first element is the maximum RE, and the second 
#' element is the maximum vector abundance and the corresponding RE
#' @export
#'
#' @examples
Calculate_max_RE_DF <- function(df, splitting_factor){
  
  split_df <- split(df, 
                    df[,splitting_factor])
  
  max_RE_DF <- do.call(rbind, lapply(split_df, function(x) x[which.max(x$RE),]))
  max_Vab_DF <- do.call(rbind, lapply(split_df, function(x) x[which.max(x$NP + x$NM),]))
                                            

 return(list(max_RE_DF, max_Vab_DF))
 
}


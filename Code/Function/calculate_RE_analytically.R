#' Calculate R_effective analytically 
#'
#' Calculate R_effective with the assumption that
#' any new infection is from either human or vector. This leads
#' to RE that is squared.
#'
#' @param param
#' @param List_x
#'
#' @return
#' @export
#'
#' @examples
Calculate_analytical_REff <- function(List_x, param) {
  
   theta_H <- param["theta_H"]
  theta_P <- param["theta_P"]
  theta_M <- param["theta_M"]
  f_P <- param["f_P"]
  f_M <- param["f_M"]
  gamma <- param["gamma"]
  mu_H <- param["mu_H"]
  c_MP <- param["c_MP"]
  c_PM <- param["c_PM"]
  c_MM <- param["c_MM"]
  c_PP <- param["c_PP"]
  ntime <- param["ntime"]

  HS <- List_x[[1]]
  HI <- List_x[[2]]
  HR <- List_x[[3]]

  PS <- List_x[[4]]
  PI <- List_x[[5]]

  MS <- List_x[[6]]
  MI <- List_x[[7]]

  NH <- HS + HI + HR
  NP <- PS + PI
  NM <- MS + MI

  Psi_P <- (c_MP * NM) + c_PP * (PS + (2 * PI))
  Psi_M <- (c_PM * NP) + c_MM * (MS + (2 * MI))
  Psi_C <- (c_MP * PI) * (c_PM * MI)
  Gamma <- (Psi_P * Psi_M) - Psi_C

  wait_time <- sqrt(((Psi_P * Psi_M) - Psi_C) * (gamma + mu_H)) * NH

  HtoH_P <- Psi_M * ((theta_H * f_P * PS) * (theta_P * f_P))
  HtoH_M <- Psi_P * ((theta_H * f_M * MS) * (theta_M * f_M))

  MtoP <- (c_PM * MI) * ((theta_H * f_P * PS) * (theta_M * f_M))
  PtoM <- (c_MP * PI) * ((theta_H * f_M * MS) * (theta_P * f_P))

  RE <- (sqrt(HS) * (sqrt(HtoH_M - PtoM + HtoH_P - MtoP))) / wait_time

  RE_DF <- cbind.data.frame(
    time = seq(1, ntime),
    RE = RE,
    NP = NP,
    NM = NM,
    PS = PS,
    MS  =MS,
    HS = HS,
    HtoM = ((theta_H * f_M) * MS) / ((gamma + mu_H) * NH),
    MtoH = ((Psi_P * (theta_M * f_M) * HS) - (c_MP * PI * (theta_P * f_P * HS))) /
      ((Gamma) * NH),
    HtoP = ((theta_H * f_P) * PS) / ((gamma + mu_H) * NH),
    PtoH = ((Psi_M * (theta_P * f_P) * HS) - (c_PM * MI * (theta_M * f_M * HS))) /
      ((Gamma) * NH),
    wait_time = wait_time,
    Gamma,
    Psi_M,
    Psi_P
  )

  return(RE_DF)
}

#' Calculate R_effective (human) analytically
#'
#'Calculate R_effective with the assuption that any new infection
#'is ONLY from the human
#'
#' @param List_x 
#' @param param 
#'
#' @return
#' @export
#'
#' @examples
Calculate_Human_REff <- function(List_x, param) {
   
   theta_H <- param["theta_H"]
   theta_P <- param["theta_P"]
   theta_M <- param["theta_M"]
   f_P <- param["f_P"]
   f_M <- param["f_M"]
   gamma <- param["gamma"]
   mu_H <- param["mu_H"]
   c_MP <- param["c_MP"]
   c_PM <- param["c_PM"]
   c_MM <- param["c_MM"]
   c_PP <- param["c_PP"]
   ntime <- param["ntime"]
   
   HS <- List_x[[1]]
   HI <- List_x[[2]]
   HR <- List_x[[3]]
   
   PS <- List_x[[4]]
   PI <- List_x[[5]]
   
   MS <- List_x[[6]]
   MI <- List_x[[7]]
   
   NH <- HS + HI + HR
   NP <- PS + PI
   NM <- MS + MI
   
   Psi_P <- (c_MP * NM) + c_PP * (PS + (2 * PI))
   Psi_M <- (c_PM * NP) + c_MM * (MS + (2 * MI))
   Psi_C <- (c_MP * PI) * (c_PM * MI)
   Gamma <- (Psi_P * Psi_M) - Psi_C
   
   wait_time <- Gamma * (gamma + mu_H) * NH^2

   Primary <-   HS * (theta_H * f_P * PS) * ((Psi_M * (theta_P * f_P) - (c_PM * MI * (theta_M * f_M))))
   Secondary <- HS*  (theta_H * f_M * MS) * ((Psi_P * (theta_M * f_M) - (c_MP * PI * (theta_P * f_P))))

   
   RE <- (Primary/wait_time) + (Secondary/wait_time) 
   
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
      wait_time = wait_time,
      Gamma,
      Psi_M,
      Psi_P
   )
   
   return(RE_DF)
}

Calculate_Human_Reff_Expanded <- function(df_expand, param){
  
  theta_H <- param["theta_H"]
  theta_P <- param["theta_P"]
  theta_M <- param["theta_M"]
  f_P <- param["f_P"]
  f_M <- param["f_M"]
  gamma <- param["gamma"]
  mu_H <- param["mu_H"]
  c_MP <- param["c_MP"]
  c_PM <- param["c_PM"]
  c_MM <- param["c_MM"]
  c_PP <- param["c_PP"]
  
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
  
  Primary <-   HS * (theta_H * f_P * PS) * ((Psi_M * (theta_P * f_P) - (c_PM * MI * (theta_M * f_M))))
  Secondary <- HS*  (theta_H * f_M * MS) * ((Psi_P * (theta_M * f_M) - (c_MP * PI * (theta_P * f_P))))
  
  
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


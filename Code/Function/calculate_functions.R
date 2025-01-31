#' Calculate the effective reproductive number (human) derived analytically
#'
#' Calculate the effective reproductive number (RE) with the assumption that any
#' new infection is ONLY from the human. The equation is derived analytically
#' (check manuscript), and we use the simulation output to track how RE changes
#' over time.
#'
#' @param List_x (list) The model output. Should have length of 7 to account
#' for the seven compartments.
#'
#' @param param (data.frame) Parameters used from simulating the model
#'
#' @return (data.frame) A that includes the RE, the total number of primary (NP)
#' and secondary vectors (NM), primary (PS) and secondary (MS) susceptibles,
#' human susceptibles (HS). The contribution to the RE from the primary vector (
#' PtoH) and the secondary vector (MtoH), and finally the 'waiting time' (wait_time)
#'
#' @export
#'
#' @examples Calculate_Human_REff(model_output, param_standard)
Calculate_Human_REff <- function(List_x, param) {
  if (length(List_x) != 7 | class(List_x) != "list") {
    stop("This is not the correct format: x`are you sure this is the model output?")
  }

  # Parameters
  theta_H <- param["theta_H"] # Transmission probability of human
  theta_P <- param["theta_P"] # Transmission probability of primary vector
  theta_M <- param["theta_M"] # Transmission probability of secondary vector
  f_P <- param["f_P"] # Primary biting rate
  f_M <- param["f_M"] # Secondary biting rate
  gamma <- param["gamma"] # Recovery rate
  mu_H <- param["mu_H"] # Human mortality rate
  mu_V <- param["mu_V"]
  c_MP <- param["c_MP"] # Competition of secondary on primary
  c_PM <- param["c_PM"] # Competition of primary on secondary
  c_MM <- param["c_MM"] # Competition of secondary on secondary
  c_PP <- param["c_PP"] # Competition of primary on primary
  ntime <- param["ntime"] # How long does the simulation run for?

  # States
  HS <- List_x[[1]] # human susceptible
  HI <- List_x[[2]] # human infected
  HR <- List_x[[3]] # human recovered

  PS <- List_x[[4]] # primary susceptible
  PI <- List_x[[5]] # primary infected

  MS <- List_x[[6]] # secondary susceptible
  MI <- List_x[[7]] # secondary infected

  NH <- HS + HI + HR # Total humans
  NP <- PS + PI # Total primary vectors
  NM <- MS + MI # Total secondary vectors

  # Related to the average dwell time
  ### Waiting time denominator
  human_wait_time <- gamma + mu_H

  Primary_HP <- (theta_H * f_P * (PS/NH)) * (1/ human_wait_time) 
  Primary_PH <- (theta_P * f_P * (HS/NH)) * (1/mu_V) 
  
  Secondary_HM <- (theta_H * f_M * (MS/NH)) * (1/ human_wait_time) 
  Secondary_MH <- (theta_M * f_M * (HS/NH)) * (1/mu_V) 
  


  RE <-   (Primary_HP * Primary_PH) + (Secondary_HM * Secondary_MH)

  RE_DF <- cbind.data.frame(
    time = seq(1, ntime),
    RE = RE,
    NP = NP,
    NM = NM,
    PS = PS,
    MS = MS,
    HS = HS,
    HI = HI,
    PtoH =   (Primary_HP * Primary_PH) ,
    MtoH = (Secondary_HM * Secondary_MH),
    Primary_HP = Primary_HP,
    Primary_PH  = Primary_PH,
    Secondary_HM = Secondary_HM,
    Secondary_MH  = Secondary_MH
  )
  return(RE_DF)
}






#' Calculate the reproductive number (human) derived analytically
#'
#' Function is different than `Calculate_Human_Reff` because instead of
#' the model output, you
#'
#' @param df_expand
#' @param param
#'
#' @returns
#' @export
#'
#' @examples
Calculate_Human_Reff_Expanded <- function(df_expand, param) {
  # Parameters
  theta_H <- param["theta_H"] # Transmission probability of human
  theta_P <- param["theta_P"] # Transmission probability of primary vector
  theta_M <- param["theta_M"] # Transmission probability of secondary vector

  f_P <- param["f_P"] # Primary biting rate
  f_M <- param["f_M"] # Secondary biting rate

  gamma <- param["gamma"] # Recovery rate

  mu_H <- param["mu_H"] # Human mortality rate

  c_MP <- param["c_MP"] # Competition of secondary on primary
  c_PM <- param["c_PM"] # Competition of primary on secondary
  c_MM <- param["c_MM"] # Competition of secondary on secondary
  c_PP <- param["c_PP"] # Competition of primary on primary

  HS <- 1000
  NH <- 1000

  NP <- df_expand["NP"]
  PI <- NP * 0.01
  PS <- NP - PI

  NM <- df_expand["NM"]
  MI <- NM * 0.01
  MS <- NM - MI


  # Related to the average dwell time
  Psi_P <- (c_MP * NM) + c_PP * (PS + (2 * PI))
  Psi_M <- (c_PM * NP) + c_MM * (MS + (2 * MI))
  Psi_C <- (c_MP * PI) * (c_PM * MI)
  Gamma <- (Psi_P * Psi_M) - Psi_C

  wait_time <- Gamma * (gamma + mu_H) * NH^2

  Primary <- HS * (theta_H * f_P * PS) * ((Psi_M * (theta_P * f_P) -
    (c_PM * MI * (theta_M * f_M))))

  Secondary <- HS * (theta_H * f_M * MS) * ((Psi_P * (theta_M * f_M) -
    (c_MP * PI * (theta_P * f_P))))


  RE <- (Primary / wait_time) + (Secondary / wait_time)

  RE_DF <- cbind.data.frame(
    RE = RE,
    NP = NP,
    NM = NM,
    PS = PS,
    MS = MS,
    HS = HS,
    PtoH = Primary / wait_time,
    MtoH = Secondary / wait_time
  )

  colnames(RE_DF) <- c("RE", "NP", "NM", "PS", "MS", "HS", "PtoH", "MtoH")

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



#' Calculate steady state (Caitlin Lienkaemper wrote this)
#'
#'Calculates when the system should reach the steady state
#'(If it doesn't work- blame her)
#'
#' @param v
#' @param max_diff
#'
#' @returns
#' @export
#'
#' @examples
Calculate_steady_state <- function(v, max_diff = 10^(-6)) {
  diffs <- diff(v)
  n <- length(diffs)
  for (i in seq_len(n)) {
    if (max(abs(diffs[i:n])) < max_diff) {
      return(c(i, v[i]))
    }
  }
  return("that aint it chief")
}

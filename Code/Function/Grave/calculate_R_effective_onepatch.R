calculate_R_effective_modified <- function(full_list, param) {
  
  theta_P <- param["theta_P"]
  theta_M <- param["theta_M"]
  theta_H <- param["theta_H"]
  gamma <- param["gamma"]
  mu_H <- param["mu_H"]
  f_P <- param["f_P"]
  f_M <- param["f_M"]
  B_V <- param["B_V"]
  c_PP <- param["c_PP"]
  c_MM <- param["c_MM"]
  c_MP <- param["c_MP"]
  c_PM <- param["c_PM"]
  ntime <- param["ntime"]

  HS <- full_list[[1]]
  HI <- full_list[[2]]
  HR <- full_list[[3]]
  PS <- full_list[[4]]
  PI <- full_list[[5]]
  MS <- full_list[[6]]
  MI <- full_list[[7]]

  time_RE <- NULL

  for (j in seq(1, ntime)) {
    NH <- HS[j] + HI[j] + HR[j]
    NP <- PS[j] + PI[j]
    NM <- MS[j] + MI[j]

    HS_ratio <- HS[j] / NH
    PS_ratio <- PS[j] / NH
    MS_ratio <- MS[j] / NH

    F_mat <- Matrix(
      c(
        0, (theta_P * f_P * HS_ratio), (theta_M * f_M * HS_ratio),
        (theta_H * f_P * PS_ratio), 0, 0,
        (theta_H * f_M * MS_ratio), 0, 0
      ),
      byrow = TRUE, ncol = 3, sparse = TRUE
    )

    F_mat[!(is.finite(F_mat))] <- 0

    V_mat <- Matrix(
      c(
        (gamma + mu_H), 0, 0,
        0, (c_MP * NM) + (c_PP * NP), 0,
        0, 0, (c_PM * NP) + (c_MM * NM)
      ),
      byrow = TRUE, ncol = 3, sparse = TRUE
    )

    FV <- (F_mat %*% solve(V_mat))
    RE <- max(eigen(FV)$values)

    time_RE[[j]] <- data.frame(
      time = j,
      RE = RE,
      N_P = NP,
      N_M = NM,
      H_S = HS[j],
      H_I = HI[j],
      P_S = PS[j],
      P_I = PI[j],
      M_I = MI[j],
      M_S = MS[j]
    )
  }

  return(do.call(rbind, time_RE))
}

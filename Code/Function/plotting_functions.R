#' Plot the total abundance of the different classes
plot_list_groups <- function(full_list) {
  compartment_labels <- c(
    "HS", "HI", "HR",
    "PS", "PI",
    "MS", "MI"
  )

  data_frame_plot <- lapply(full_list, function(x) {
    data.frame(x, time = seq_len(nrow(x)))
  })

  # Assigning the classes so we can combine into one big data.frame
  for (i in seq(1, length(data_frame_plot))) {
    data_frame_plot[[i]]$id <- compartment_labels[i]
  }

  final_data_frame <- do.call(rbind, data_frame_plot)

  final_melted_frame <- melt(
    final_data_frame,
    id.vars = list("time", "id")
  )

  final_melted_frame$id <- factor(
    final_melted_frame$id,
    levels = compartment_labels
  )

  full_plot_GG <- ggplot(
    final_melted_frame,
    aes(x = time, y = value, color = variable)
  ) +
    geom_line() +
    facet_wrap(~id, scales = "free_y") +
    xlab("Time") +
    ylab("Abundance") +
    theme_classic() +
    scale_color_viridis(discrete = TRUE, option = "turbo") +
    theme(
      legend.position = "none",
      strip.background = element_blank()
    )

  return(full_plot_GG)
}



# 02-Figure
#' Plot the RE over time as well as side panels,
plot_RE_dynamics <- function(df, dstb, postdisturb = "No", RE_limit) {
  if (any((colnames(df) %in% c("NP", "NM"))) == FALSE) {
    stop("You are missing NP and/or NM, check column names")
  }

  if (any((colnames(df) %in% c("RE"))) == FALSE) {
    stop("No RE, check column names ")
  }

  if (postdisturb == "No") {
    subset_df <- df
    subset_df$time <- subset_df$time - dstb
    maximum_RE <- Calculate_max_RE_DF(subset_df, "id")[[1]]
    minimum_RE <- Calculate_max_RE_DF(subset_df, "id")[[2]]
  } else if (postdisturb == "Yes") {
    subset_df <- df
    maximum_RE <- Calculate_max_RE_DF(subset_df, "id")[[1]]
  }

  time_RE_GG <-
    ggplot(subset_df, aes(x = time)) +
    geom_path(aes(
      y = (RE),
      group = as.factor(1 - id),
      color = as.factor(1 - id)
    ), linewidth = 0.3) +
    geom_hline(
      data = subset(subset_df, subset_df$time == 0),
      aes(yintercept = (RE)), color = "black", alpha = 0.7, linetype = 3
    ) +
    scale_color_discrete_sequential(palette = "ag_Sunset", rev = FALSE) +
    xlab("Time since disturbance") +
    ylab(expression("Reproductive number (" * R[t] * ")")) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      aspect.ratio = 1
    )


  return(time_RE_GG)
}
#' Plot the total vector abundance over time with the RE being the color.
plot_NV_RE <- function(df, dstb, postdisturb = "No", RE_limit) {
  if (any((colnames(df) %in% c("NP", "NM"))) == FALSE) {
    stop("You are missing NP and/or NM, check column names")
  }

  if (any((colnames(df) %in% c("RE"))) == FALSE) {
    stop("No RE, check column names ")
  }

  if (postdisturb == "No") {
    subset_df <- df
    subset_df$time <- subset_df$time - dstb
    maximum_RE <- Calculate_max_RE_DF(subset_df, c("param", "id"))[[1]]
  } else if (postdisturb == "Yes") {
    subset_df <- df
    maximum_RE <- Calculate_max_RE_DF(subset_df, c("param", "id"))[[1]]
  }


  NV_RE_GG <-
    ggplot(subset_df, aes(x = time)) +
    geom_path(aes(y = (NP + NM), color = RE, group = id), linewidth = 0.3) +
    geom_point(
      data = maximum_RE,
      aes(x = time, y = (NP + NM), group = id, color = RE),
      size = 1,
    ) +
    geom_point(
      data = maximum_RE,
      aes(x = time, y = (NP + NM), group = id),
      shape = 21, size = 1.2, fill = NA
    ) +
    geom_hline(
      data = subset(subset_df, subset_df$time == 0),
      aes(yintercept = NP + NM), color = "black", alpha = 0.7, linetype = 3
    ) +
    xlab("Time since disturbance") +
    ylab("Total vector abundance") +
    scale_color_viridis(name = expression(R[t])) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black")
    )
  return(NV_RE_GG)
}
#' Plot trajectories of primary and secondary vectors with the RE
#' as color
plot_NP_NM_RE <- function(df, postdisturb, RE_limit) {
  if (any((colnames(df) %in% c("NP", "NM"))) == FALSE) {
    stop("You are missing NP and/or NM, check column names")
  }
  if (any((colnames(df) %in% c("RE"))) == FALSE) {
    stop("No RE, check column names ")
  }

  if (postdisturb == "No") {
    subset_df <- df
    subset_df$time <- subset_df$time - 9125
    maximum_RE <- Calculate_max_RE_DF(subset_df, c("param", "id"))[[1]]
  } else if (postdisturb == "Yes") {
    subset_df <- df
    maximum_RE <- Calculate_max_RE_DF(subset_df, c("param", "id"))[[1]]
  }

  NP_NM_RE_GG <- ggplot(
    subset_df, aes(x = NP, y = NM)
  ) +
    geom_path(
      aes(color = RE),
      linewidth = 0.3,
      linejoin = "round",
      lineend = "round"
    ) +
    geom_point(
      data = subset(
        subset_df,
        subset_df$time == 0
      ),
      aes(x = NP, y = NM), color = "black"
    ) +
    geom_point(
      data = maximum_RE,
      aes(x = NP, y = NM, group = id, color = RE),
      size = 1
    ) +
    geom_point(
      data = maximum_RE,
      aes(x = NP, y = NM, group = id),
      shape = 21, size = 1.2, fill = NA
    ) +
    scale_color_viridis(name = expression(R[t])) +
    xlab(expression("# of primary vectors " * "(" * N^P * ")")) +
    ylab(expression("# of secondary vectors " * "(" * N^M * ")")) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black")
    )

  return(NP_NM_RE_GG)
}



# 03-Figure
#' Plot the maximum secondary contribution of RE
plot_NM_REff <- function(df, postdisturb) {
  if (any((colnames(df) %in% c("NP", "NM"))) == FALSE) {
    stop("You are missing NP and/or NM, check column names")
  }

  if (any((colnames(df) %in% c("RE"))) == FALSE) {
    stop("No RE, check column names ")
  }

  if (postdisturb == "No") {
    subset_df <- df
    subset_df$time <- subset_df$time - 9125
    maximum_RE <- Calculate_max_RE_DF(subset_df, "id")[[1]]
  } else if (postdisturb == "Yes") {
    subset_df <- df
    maximum_RE <- Calculate_max_RE_DF(subset_df, "id")[[1]]
  }


  NM_Reff_GG <- ggplot(
    subset_df, aes(x = time, y = MtoH / RE, group = id, color = RE)
  ) +
    geom_line(
      linewidth = 0.9,
      linejoin = "round",
      lineend = "round"
    ) +
    geom_point(
      data = maximum_RE,
      aes(x = time, y = MtoH / RE, group = id)
    ) +
    xlab("Time since disturbance") +
    ylab("Secondary vector contribution to Rt") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 15, color = "black")
    )

  return(NM_Reff_GG)
}

plot_m_contribution_heatmap <- function(df) {
  param_standard <- get_parameters("standard")
  theta_M_standard <- param_standard["theta_M"]
  f_P_standard <- param_standard["f_P"]
  theta_P_standard <- param_standard["theta_P"]



  GG_heat_mtoH <- ggplot(
    df,
    aes(
      x = (f_M * theta_M_standard) / (f_P_standard * theta_P_standard),
      y = 1 - mortality_P, z = max_mtoH
    )
  ) +
    geom_contour_filled(color = "black", linewidth = 0.7, binwidth = 0.10) +
    geom_textcontour(bins = 10, size = 2.5, textcolour = "white") +
    coord_equal() +
    xlab("Transmission efficiency of secondary vector \n in relation to primary vector") +
    ylab("Proportion of primary vectors removed") +
    scale_fill_viridis(option = "mako", name = expression("Secondary vector \ncontribution to" ~ R[t]), discrete = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    ) +
    guides(color = guide_legend(override.aes = list(size = 0.5)))




  GG_heat_RE <- ggplot(
    df,
    aes(
      x = (f_M * theta_M_standard) / (f_P_standard * theta_P_standard),
      y = 1 - mortality_P, z = max_RE
    )
  ) +
    geom_contour_filled(color = "black", linewidth = 0.7, binwidth = 5) +
    coord_equal() +
    xlab("Transmission efficiency of secondary vector \n in relation to primary vector") +
    ylab("Proportion of primary vectors removed") +
    scale_fill_viridis(option = "mako", name = expression("Increase in" ~ R[t]), discrete = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    theme(
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )


  return(list(GG_heat_mtoH, GG_heat_RE))
}
#' Plot the secondary contribution of RE overtime time
plot_m_contribution_lineplot <- function(df, dstb_time) {
  GG_line_RE <- ggplot(
    df,
    aes(
      x = time - dstb_time,
      y = max_mtoH,
      linetype = as.factor(standardized_ratio),
      label = as.factor(standardized_ratio)
    )
  ) +
    geom_line() +
    geom_point(
      data = subset(
        df,
        df$time == dstb_time - 1
      ),
      aes(x = time - dstb_time, y = max_mtoH)
    ) +
    scale_color_grey() +
    geom_dl(
      method = list(dl.combine("last.points")), cex = 0.1
    ) +
    xlab("Time since disturbance") +
    ylab("Secondary vector \ncontribution to" ~ R[t]) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      legend.position = "none",
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, color = "black")
    )


  return(GG_line_RE)
}





#' Plot primary vectors, secondary vectors, and human infections
plot_state_dynamics <- function(df, mortality_P = 0.25) {
  if (any((colnames(df) %in% c("NP", "NM"))) == FALSE) {
    stop("You are missing NP and/or NM, check column names")
  }

  if (any((colnames(df) %in% c("RE"))) == FALSE) {
    stop("No RE, check column names ")
  }

  subset_df <- subset(df, df$id == mortality_P & df$time < 12000)

  state_GG_V <-
    ggplot(subset_df, aes(x = time - 9125)) +
    geom_line(aes(y = NP, color = "NP"), linewidth = 1.1) +
    geom_line(aes(y = NM, color = "NM"), linewidth = 1.2) +
    scale_color_manual(
      name = "Vector",
      breaks = c("NP", "NM"),
      label = c("Primary", "Secondary"),
      values = c("NP" = "#646198", "NM" = "#D65739")
    ) +
    xlab("Time") +
    ylab("Abundance") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 15, color = "black"),
      legend.position = "top"
    )

  state_GG_HI <-
    ggplot(
      subset_df,
      aes(x = time - 9125, y = NewHumanCases, group = 1)
    ) +
    geom_line(color = "black", linewidth = 1.1) +
    geom_line(
      data = subset_df,
      aes(x = time - 9125, y = NewHumanCases), group = 1
    ) +
    xlab("Time") +
    ylab("New human \n infections") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 15, color = "black")
    )


  return(state_GG_V / state_GG_HI)
}

#' Plot the total vector abundance for disturbance intensity
plot_comparison_RE <- function(df, split_variable) {
  if (any((colnames(df) %in% c("NP", "NM"))) == FALSE) {
    stop("You are missing NP and/or NM, check column names")
  }

  if (any((colnames(df) %in% c("RE"))) == FALSE) {
    stop("No RE, check column names ")
  }

  max_RE_ab_list <- Calculate_max_RE_DF(df, split_variable)

  max_RE_ab_GG <- ggplot(
    max_RE_ab_list[[1]],
    aes(x = as.factor(1 - id), y = NP + NM, fill = RE)
  ) +
    geom_point(
      shape = 21,
      size = 3
    ) +
    geom_point(
      data = max_RE_ab_list[[3]],
      aes(
        x = as.factor(1 - id),
        y = NP + NM,
        fill = RE
      ),
      shape = 22,
      size = 3
    ) +
    xlab("Disturbance intensity on primary vector") +
    ylab("Total vector abundance") +
    scale_fill_continuous_divergingx(
      name = expression(R[t]),
      mid = 1, n_interp = 11, palette = "Roma",
      rev = TRUE, limits = c(0.8, 1.6)
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 15, color = "black")
    )


  return(max_RE_ab_GG)
}



#' Plot heatmap with different intensity
plot_heatmapintensity <- function(df) {
  ggplot(
    Maximum_RE_Mort_PM,
    aes(
      x = as.factor(Mort_P),
      y = as.factor(Mort_M),
      fill = RE
    )
  ) +
    geom_raster() +
    scale_fill_viridis(option = "viridis", name = expression(R[t])) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    xlab(expression("Disturbance intensity of primary vector " * "(" * mu[P] * ")")) +
    ylab(expression("Disturbance intensity of secondary vector " * "(" * mu[M] * ")")) +
    coord_equal() +
    theme(
      axis.text = element_text(color = "black", size = 13),
      axis.title = element_text(color = "black", size = 14)
    )
}

plot_competition_PM <- function(df) {
  Panel_1_GG <-
    ggplot(df, aes(
      x = standardized,
      y = (RE),
      group = as.factor(1 - mortality_P)
    )) +
    geom_line(aes(color = as.factor(1 - mortality_P))) +
    geom_rect(
      xmin = 1.5,
      xmax = 2,
      ymin = 0,
      ymax = 3.3,
      fill = "grey",
      alpha = 0.01
    ) +
    geom_rect(
      xmin = 0.5,
      xmax = .6766667,
      ymin = 0,
      ymax = 3.3,
      fill = "grey",
      alpha = 0.01
    ) +
    geom_point(
      data = subset(df, df$standardized == 1.00),
      aes(
        x = standardized, y = RE,
        group = as.factor(1 - mortality_P)
      ),
      size = 1.5
    ) +
    geom_point(
      data = subset(df, df$standardized %in% c(0.75, 1.25) &
        df$mortality_P == 0.01),
      aes(x = standardized, y = RE), size = 1.5, color = "black"
    ) +
    scale_color_discrete_sequential(
      palette = "ag_Sunset", rev = FALSE,
      n = 5
    ) +
    geom_vline(xintercept = c(.6766667)) +
    geom_vline(xintercept = c(1.5)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, 2)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(expression("Multiplier of primary on secondary competition (" * c[PM] * ")")) +
    ylab(expression("Increase from baseline " * R[t]^"*")) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black")
    )

  Panel_2_GG <- ggplot(
    df, aes(
      x = c_PM / c_PM_standard,
      y = (max_NM),
      group = as.factor(1 - mortality_P)
    )
  ) +
    geom_line(aes(color = as.factor(1 - mortality_P))) +
    geom_rect(
      xmin = 1.5,
      xmax = 2,
      ymin = 0,
      ymax = 1250,
      fill = "grey",
      alpha = 0.01
    ) +
    geom_rect(
      xmin = 0.5,
      xmax = .6766667,
      ymin = 0,
      ymax = 1250,
      fill = "grey",
      alpha = 0.01
    ) +
    geom_point(
      data = subset(df, df$standardized == 1),
      aes(x = standardized, y = max_NM), size = 1.5, color = "black"
    ) +
    geom_point(
      data = subset(df, df$standardized %in% c(0.75, 1.25) &
        df$mortality_P == 0.01),
      aes(x = standardized, y = max_NM), size = 1.5, color = "black"
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, 2)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_discrete_sequential(palette = "ag_Sunset", rev = FALSE, n = 5) +
    geom_vline(xintercept = c(1.5)) +
    geom_vline(xintercept = c(.6766667)) +
    xlab(expression("Multiplier of primary on secondary competition (" * c[PM] * ")")) +
    ylab(expression("Increase from " * N[M]^"*")) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black")
    )

  return(list(Panel_1_GG, Panel_2_GG))
}

plot_phaseplot <- function(df_x, isocline_df) {
  GG_phaseplot <- ggplot() +
    geom_segment(data = isocline_df$Isocline[[1]], aes(x = x, xend = xend, y = y, yend = yend), color = "#646198") +
    geom_segment(data = isocline_df$Isocline[[2]], aes(x = x, xend = xend, y = y, yend = yend), color = "#D65739") +
    geom_point(aes(x = df_x[1, "NP"], y = df_x[1, "NM"])) +
    theme_classic() +
    scale_color_viridis(name = "Rt") +
    geom_quiver(
      data = isocline_df$Phaseplot,
      aes(
        x = Var1,
        y = Var2,
        u = slope_P,
        v = slope_M, color = RE
      ),
      vecsize = 2,
      alpha = 1
    ) +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 2000)) +
    geom_path(data = df_x, aes(x = NP, y = NM), size = 0.5, alpha = 1) +
    xlab("Primary vectors") +
    ylab("Secondary vectors") +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.4, "cm")
    )

  return(GG_phaseplot)
}

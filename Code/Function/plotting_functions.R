#' Plot the total abundance of the different classes
#'
#' When given the results of `discrete_trito_model`, 
#' plot out a ggplot with 7 facets (one for each class).
#' On the x-axis is time and the abundance on the y-axis.
#'
#' @param full_list (list) The full model output that comes from 
#' direct simulation - not processed (so no RE)
#'
#' @return A ggplot object
#' @export
#' @examples
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
    aes(x = time, y = (value), color = variable)
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

#' Plot trajectories of primary and secondary vectors with the RE
#' as color
#'
#' After a disturbance,plot out how the trajectory of the total primary
#' and secondary vector with the RE as color. On the x-axis is the total
#' primary vector, on the y-axis is the total secondary vector. The 
#' different lines represent the different disturbance intensity 
#' 
#'
#' @param df (data.frame)
#'
#' @returns
#' @export
#'
#' @examples
plot_NP_NM_RE <- function(df) {
  if (any((colnames(df) %in% c("NP", "NM"))) == FALSE) {
    stop("You are missing NP and/or NM, check column names")
  }

  if (any((colnames(df) %in% c("RE"))) == FALSE) {
    stop("No RE, check column names ")
  }

  NP_NM_RE_GG <- ggplot(
    df, aes(x = NP, y = NM, color = RE)
  ) +
    geom_path(
      linewidth = 1.2,
      linejoin = "round",
      lineend = "round"
    ) +
    geom_point(
      data = subset(
        RE_mortality_P_post,
        RE_mortality_P_post$time == 9125
      ),
      aes(x = NP, y = NM), color = "black"
    ) +
    scale_colour_gradient2(
      low = "#448C81",
      mid = "#EFEAD7",
      high = "#C25D31",
      midpoint = 1,
      name = expression(R[E])
    ) +
    xlab(expression("Abundance of primary vectors " * "(" * N[P] * ")")) +
    ylab(expression("Abundance of secondary vectors " * "(" * N[M] * ")")) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 15, color = "black")
    )
  return(NP_NM_RE_GG)
}

#' Plot the total abundance of the different classes
#'
#' When given the results of `discrete_trito_model`, this
#' will plot out 7 facets (one for each class) with
#' time on the x-axis and the abundance on the y-axis.
#'
#' @param full_list The full model output that comes from the
#'
#' @return A ggplot object
#' @export
#'
#' @examples
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
    aes(x = time, y = (value), color = variable)
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

#' Plot the total vector abundance over time with the RE being the color.
#'
#' Shows how the total vector abundance (both primary and secondary) changes
#' over time as well as the RE. Points indicate the maximum RE reached and is
#' used to show how the maximum vector abundance does not correspond to the maximum
#' RE. The x-axis is time, the y-axis is the total vector abundance.

#' @param df 
#'
#' @returns A ggplot item
#' @export
#'
#' @examples
plot_NV_RE <- function(df, postdisturb = "No") {
        
  if (any((colnames(df) %in% c("NP", "NM"))) == FALSE) {
    stop("You are missing NP and/or NM, check column names")
  }

  if (any((colnames(df) %in% c("RE"))) == FALSE) {
    stop("No RE, check column names ")
  }

  if (postdisturb == "No"){
  subset_df <- subset(df, df$time > 9370)  
  maximum_RE <- Calculate_max_RE_DF(subset_df , "id")[[1]]
  } else if (postdisturb == "Yes"){
  subset_df <- df
  maximum_RE <- Calculate_max_RE_DF(subset_df , "id")[[1]]    
  }
  
  
  NV_RE_GG <-
    ggplot(subset_df, aes(x = time )) +
    geom_line(aes(y = (NP + NM), color = RE, group = id), linewidth = 0.9) +
    geom_point(data = maximum_RE, 
      aes(x = time , y = (NP+ NM), group = id), shape = 5) + 
    xlab("Time since disturbance") +
    ylab("Total vector abundance") +
    scale_colour_gradient2(
      low = "#448C81",
      mid = "#EFEAD7",
      high = "#C25D31",
      midpoint = 1,
      name = expression(R[E])) + 
    theme_classic() +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 15, color = "black"),
      legend.position = "top"
    )
  return(NV_RE_GG)
}


#' Plot primary vectors, secondary vectors, and human infections
#'
#' Plot the total abundance of both the primary and secondary 
#' vector 
#'
#' @param df
#'
#' @returns
#' @export
#'
#' @examples
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
#'
#' Plot the total vector abundance for when the RE is maximized for 
#' different disturbance intensity (circles) AND the correpsonding 
#' RE when the total vector abundance is maximized. On the 
#' x-axis we have the disturbance intensity, the y-axis as the 
#' total vector abundance (NP + NM)
#'
#'
#' @param df
#'
#' @returns
#' @export
#'
#' @examples
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
    aes(x = 1 - id, y = NP + NM, fill = RE)
  ) +
    geom_point(
      shape = 21,
      size = 3
    ) +
    geom_point(
      data = max_RE_ab_list[[2]],
      aes(
        x = 1 - id,
        y = NP + NM,
        fill = RE
      ),
      shape = 22,
      size = 3
    ) +
    xlab("Mortality of primary vector") +
    ylab("Total vector abundance") +
    scale_fill_viridis() +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15)
    )


  return(max_RE_ab_GG)
}


#' Plot the heat-map of the R0 depending on the abundance of the
#' primary and secondary vectors.
#'
#' @param df The data.frame with the total primary vectors
#' and the secondary vectors.
#'
#' @returns A heatmap plot is the
#' @export
#'
#' @examples
plot_heatmapR0 <- function(df) {
  if (any((colnames(df) %in% c("NP", "NM"))) == FALSE) {
    stop("You are missing NP and/or NM, check column names")
  }

  if (any((colnames(df) %in% c("RE"))) == FALSE) {
    stop("No RE, check column names ")
  }

  heatmap_GG <- ggplot(
    df,
    aes(x = NP, y = NM, fill = RE, color = RE)
  ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_tile() +
    xlab(expression("Abundance of primary vectors " * "(" * N[P] * ")")) +
    ylab(expression("Abundance of secondary vectors " * "(" * N[M] * ")")) +
    scale_fill_viridis(
      option = "rocket",
      name = expression(R[E])
    ) +
     scale_color_viridis(
      option = "rocket",
      name = expression(R[E])
    ) +
    theme(
      axis.text = element_text(size = 12.5, color = "black"),
      axis.title = element_text(size = 14)
    )


  return(heatmap_GG)
}

#' Plot heatmap with different intensity
#'
#' @param df
#'
#' @returns
#' @export
#'
#' @examples
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
    scale_fill_viridis(option = "rocket", name = expression(R[E])) +
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

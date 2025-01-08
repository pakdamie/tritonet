# Interspecific competition

c_PP_standard <- 4.5e-06 ## Competition effect of p.vectors
c_MM_standard <- 3e-06 ## Competition effect of s.vectors
modifier <- seq(0.5, 1.5, by = 0.1)


c_PP <- c_PP_standard * modifier
c_MM <- c_MM_standard * modifier
Mortality_P <- c(0.01)

intra_competition <- data.frame(expand.grid(
  c_PP, c_MM,
  Mortality_P
))
colnames(intra_competition) <- c("c_PP", "c_MM", "mortality_P")


# stupid way that i'm gonna fix

param_standard_list <- NULL
for (i in seq(1:nrow(intra_competition))) {
  param_copy <- param_standard

  param_copy["c_PP"] <- intra_competition[i, "c_PP"]
  param_copy["c_MM"] <- intra_competition[i, "c_MM"]
  param_copy["mortality_P"] <- intra_competition[i, "mortality_P"]

  param_standard_list[[i]] <- param_copy
}


RE_CM <-
  Simulate_Model_Output(
    param_standard, c("c_PP", "c_MM", "mortality_P"),
    intra_competition
  )



RE_DF <- NULL
for (k in 1:length(RE_CM)) {
  RE_tmp <- Calculate_Human_REff(RE_CM[[k]], param_standard_list[[k]])
  post_df <- na.omit(subset(RE_tmp, RE_tmp$time > 9124))

  if (nrow(post_df) == 0) {
    RE_max <- cbind(intra_competition[k, ], RE = NA, MtoH = NA)
  } else {
    RE_max <- cbind(intra_competition[k, ], RE = max(post_df$RE), 
                    MtoH = post_df[which.max(post_df$MtoH),]$MtoH/
                            post_df[which.max(post_df$MtoH),]$RE)
  }
  RE_DF[[k]] <- RE_max
}

RE_DF <- do.call(rbind, RE_DF)

intra_GG <- ggplot(
  RE_DF,
  aes(
    x = as.factor((c_PP / c_PP_standard)),
    y = as.factor((c_MM / c_MM_standard)),
    fill = RE
  )
) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("B. Intraspecific competition") +
  scale_fill_continuous_divergingx(
    name = expression(R[E]),
    mid = 1, n_interp = 11,
    palette = "Roma", rev = TRUE, limits = c(0.95, 5.25)
  ) +
  xlab(expression("Modifier of primary on primary competition" ~ (italic(c[PP])))) +
  ylab(expression("Modifier of secondary on secondary competition" ~ (italic(c[MM])))) +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black")
  )



mtoh_intra_GG <- ggplot(
  RE_DF,
  aes(
    x = as.factor((c_PP / c_PP_standard)),
    y = as.factor((c_MM / c_MM_standard)),
    fill = MtoH,
  )
) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("B. Intraspecific competition") +
  scale_fill_viridis() + 
  xlab(expression("Modifier of primary on primary competition" ~ (italic(c[PP])))) +
  ylab(expression("Modifier of secondary on secondary competition" ~ (italic(c[MM])))) +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black")
  )

                                 


(inter_GG + intra_GG) /
(mtoH_inter_GG  + mtoh_intra_GG) + 
coord_equal() +  plot_layout(guides = 'collect')


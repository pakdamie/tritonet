# Interspecific competition
c_PM_standard <- 4e-6 ## Competition effect of p.vector on s.vector
c_MP_standard <- 2e-6 ## Competition effect of s.vector on p.vector

modifier <- seq(0.5, 1.5, by = 0.1)
c_PM <- c_PM_standard * modifier
c_MP <- c_MP_standard * modifier
Mortality_P <- c(0.01)

inter_competition <- data.frame(expand.grid(c_PM, c_MP, Mortality_P))
colnames(inter_competition) <- c("c_PM", "c_MP", "mortality_P")

# stupid way that i'm gonna fix

param_standard_list <- NULL
for (i in seq(1:nrow(inter_competition))) {
  param_copy <- get_parameters("standard")

  param_copy["c_PM"] <- inter_competition[i, "c_PM"]
  param_copy["c_MP"] <- inter_competition[i, "c_MP"]
  param_copy["mortality_P"] <- inter_competition[i, "mortality_P"]

  param_standard_list[[i]] <- param_copy
}


RE_CM <-
  Simulate_Model_Output(
    get_parameters("standard"), c("c_PM", "c_MP", "mortality_P"),
    inter_competition
  )



RE_DF_inter <- NULL
for (k in 1:length(RE_CM)) {
  RE_tmp <- Calculate_Human_REff(RE_CM[[k]], param_standard_list[[k]])
  post_df <- na.omit(subset(RE_tmp, RE_tmp$time > 9124))

  if (nrow(post_df) == 0) {
    RE_max <- cbind(inter_competition[k, ], RE = NA, MtoH = NA)
  } else {
    RE_max <- cbind(inter_competition[k, ], RE = max(post_df$RE),  
                   MtoH = post_df[which.max(post_df$MtoH),]$MtoH/
                          post_df[which.max(post_df$MtoH),]$RE)
  }
  RE_DF_inter[[k]] <- RE_max
}

RE_DF_inter <- do.call(rbind, RE_DF_inter)

inter_GG <- ggplot(
  RE_DF_inter,
  aes(
    x = as.factor(c_PM / c_PM_standard),
    y = as.factor(c_MP / c_MP_standard),
    fill = RE
  )
) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("A. Interspecific competition") + 
  scale_fill_continuous_divergingx(
    name = expression(R[E]),
    mid = 1, n_interp = 11,
    palette = "Roma", rev = TRUE, limits = c(0.95, 5.25)
  ) +
  xlab(expression("Modifier of primary on secondary competition"~(italic(c[PM])))) +
  ylab(expression("Modifier of secondary on primary competition"~(italic(c[MP])))) + 
  theme(axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 14, color = 'black')) 



mtoH_inter_GG <- ggplot(
  RE_DF_inter,
  aes(
    x = as.factor(c_PM / c_PM_standard),
    y = as.factor(c_MP / c_MP_standard),
    fill = MtoH
  )
) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("A. Interspecific competition") + 
  scale_fill_viridis() +
  xlab(expression("Modifier of primary on secondary competition"~(italic(c[PM])))) +
  ylab(expression("Modifier of secondary on primary competition"~(italic(c[MP])))) + 
  theme(axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 14, color = 'black')) 

                                 
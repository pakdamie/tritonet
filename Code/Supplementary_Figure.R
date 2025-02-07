# Supplementary Codes:
RE_mortality_P_post <- readRDS("Output/df_RE_mortality_P_R0.rds")

split_Mort_P <- split(RE_mortality_P_post, RE_mortality_P_post$id)

mort_P_PM <- do.call(rbind.data.frame, lapply(
  split_Mort_P,
  function(x) {
    cbind(
      mortality_P = unique(x$id),
      NM = max(x$NM) - x[1, ]$NM,
      NP = max(x$NP) - x[1, ]$NP
    )
  }
)) |>
  melt(, id.vars = "mortality_P")


MS_abundance_GG <-
  ggplot(mort_P_PM, aes(x = 1 - mortality_P, y = value, color = variable)) +
  geom_point(size = 2.5) +
  geom_line() +
  xlab("Proportion of primary vector removed") +
  ylab("Change from equilibrium abundance") +
  scale_color_manual(
    name = "Species",
    values = c("NM" = "#646198", "NP" = "#D65739"),
    labels = c("Secondary", "Primary")
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black")
  )

ggsave(here("Figures", "MS_abundance_GG.pdf"), width = 6.5, height = 5, units = "in")

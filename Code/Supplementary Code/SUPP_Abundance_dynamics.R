# SUPP

# Retrieve "standard" parameters and set the different values of disturbances
param_standard <- get_parameters("standard")
Mortality_P <- data.frame(Mort_P = c(0.01, 0.5, 0.75))

# Simulate model output and calculate the RE
RE_mortality_P <-
  Simulate_Model_Output(param_standard,
    infection_start = "No",
    "mortality_P",
    Mortality_P
  ) |>
  lapply(Calculate_Human_Reff_Expanded, param = param_standard)

# Assign mortality_P value to each list element
for (i in seq_along(RE_mortality_P)) {
  RE_mortality_P[[i]]$id <- Mortality_P[i, ]
}
RE_mortality_P_DF <- do.call(rbind, RE_mortality_P)


RE_mortality_P_sub <- RE_mortality_P_DF[,c("time","NP","NM","id")] |>
        melt(id.vars = c("time","id"))


ggplot(RE_mortality_P_sub , aes(x = time, y = log10(value),color = variable)) +
  geom_line(size = 0.9) +
  scale_color_manual(name = "Species",
          values = c("NP" = "#646198",
                                "NM" =  "#D65739"),
                     labels = c("NP" = "Primary",
                                "NM" = "Secondary"))+
  facet_wrap(~(1-id),) +  
  theme_classic() + 
  xlab("Time") + 
  ylab("Vector abundance (log10)")+
  theme(axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 15, color = 'black'),
        strip.text = element_text(size = 16, color = 'black'),
        strip.background = element_blank())

ggsave(here("Figures", "SUPP_abundance_dynamics.pdf"),
       width = 12, height = 5, units = "in"
)

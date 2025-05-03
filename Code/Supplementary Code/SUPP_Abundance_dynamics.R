# -------------------------------------------------------------
# What does the vector abundance look like pre-disturbance,
# disturbance, and post-disturbance
# -------------------------------------------------------------
# Define the different parameter scenarios for simulation
param_standard <- get_parameters("standard")

# Define different mortality rates for the primary vector
Mortality_P <- data.frame(Mort_P = c(0, 0.25, 0.5, 0.75))

# When is the system disturbed?
dstb_time <- get_parameters("standard")[["disturbance_time"]]

# --------------------------------------------------------------------------------
# Step 1: Simulate Model Output and Calculate Re
# --------------------------------------------------------------------------------


# Iterate over each parameter set and simulate model output.
# Then, compute the expanded human effective reproductive number (Re).
RE_mortality_P <- Simulate_Model_Output(param_standard,
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


# --------------------------------------------------------------------------------
# Step 2: Subset Time Period of Interest
# --------------------------------------------------------------------------------

# Extract a specific time window from the simulation results.
# The subset includes data slightly before and after the disturbance event.

# Get the time, total primary and secondary, different disturbance
# intensity and for melting
RE_mortality_P_sub <- RE_mortality_P_DF[, c("time", "NP", "NM", "id")] |>
  melt(id.vars = c("time", "id"))

RE_mortality_P_sub <- subset(
  RE_mortality_P_sub,
  RE_mortality_P_sub$time > dstb_time - 50 &
    RE_mortality_P_sub$time < dstb_time + 300
)



# -------------------------------------------------------------
# Step 3: Plot
# -------------------------------------------------------------

ggplot(RE_mortality_P_sub,
  aes(
    x = time - dstb_time,
    y = log10(value + 1),
    color = variable
  )
) +
  geom_line(size = 0.5) +
  scale_color_manual(
    name = "Species",
    values = c(
      "NP" = "#646198",
      "NM" = "#D65739"
    ),
    labels = c(
      "NP" = "Primary",
      "NM" = "Secondary"
    )
  ) +
  facet_wrap(~ (1 - id)) +
  theme_classic() +
  xlab("Time since control") +
  ylab("Vector abundance (log10)") +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", 
                                fill = NA)
    
  )

ggsave(here("Main_Figures", "SUPP_abundance_dynamics.pdf"),
  width = 6, height = 5, units = "in"
)

# -------------------------------------------------------------
# Does perturbation only lead to the secondary vector increasing
# from its equilibrium?
# -------------------------------------------------------------


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

ggsave(here("Main_Figures", "SUPP_over_Eq_M.pdf"),
       width = 6, height = 5, units = "in"
)


# -------------------------------------------------------------
# How do varying levels of control of the primary and
# secondary vector
# -------------------------------------------------------------

# Define the different parameter scenarios for simulation
param_standard <- get_parameters("standard")

# When is the system disturbed?
dstb_time <- get_parameters("standard")[["disturbance_time"]]


Mortality_P <- c(seq(0, 1, 0.1))
Mortality_M <- c(seq(0, 1, 0.1))
Mortality_PM_Grid <- expand.grid(Mortality_P, Mortality_M)

# --------------------------------------------------------------------------------
# Step 1: Simulate Model Output and Calculate Re
# --------------------------------------------------------------------------------


# Iterate over each parameter set and simulate model output.
# Then, compute the expanded human effective reproductive number (Re).

Mod_Mort_PM <- Simulate_Model_Output(
  param_standard,
  infection_start = "No", c("mortality_P", "mortality_M"),
  Mortality_PM_Grid
)

### Calculate human RE
RE_Mort_PM <- lapply(Mod_Mort_PM, Calculate_Human_Reff_Expanded, param = param_standard)

# Assign mortality_P value to each list element

for (i in 1:nrow(Mortality_PM_Grid)) {
  RE_Mort_PM[[i]]$Mort_P <- Mortality_PM_Grid[i, 1]
  RE_Mort_PM[[i]]$Mort_M <- Mortality_PM_Grid[i, 2]
}
RE_Mort_PM_DF <- do.call(rbind, RE_Mort_PM)

### Subset time period of interest
RE_Mort_PM_DF <- subset(
  RE_Mort_PM_DF,
  RE_Mort_PM_DF$time > dstb_time - 100 &
    RE_Mort_PM_DF$time < dstb_time + 1000
)


### What is the RE above the baseline?
Maximum_RE_Mort_PM <- do.call(
  rbind,
  by(RE_Mort_PM_DF,
    RE_Mort_PM_DF[, c("Mort_P", "Mort_M")],
    function(x) {
      cbind.data.frame(
        Mort_P = x$Mort_P,
        Mort_M = x$Mort_M,
        RE = x[which.max(x$RE), ]$RE - x[x$time == dstb_time - 1, ]$RE
      )
    },
    simplify = TRUE
  )
)

# -------------------------------------------------------------
# Step 2: Plot
# -------------------------------------------------------------


ggplot(
  Maximum_RE_Mort_PM,
  aes(x = (1 - Mort_P), y = (1 - Mort_M), fill = RE)
) +
  geom_tile() +
  xlab("Proportion of primary vector removed") +
  ylab("Proportion of secondary vector removed") +
  scale_fill_viridis(expression("Increase above " * R[E]^"*")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_equal() +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black")
  )

ggsave("Figures/SUPP_Intensity_PM.pdf", units = "in", height = 6, width = 6)

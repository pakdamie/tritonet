# This is the pipeline for the entirety of the manuscript.
# Maintained by Damie Pak
library(here)
source(here("Code", "Function", "packages.R"))
sourceCpp(here("Code", "Function", "model_vectors_host.cpp"))
source(here("Code", "Function", "simulate_functions.R"))
source(here("Code", "Function", "calculate_functions.R"))
source(here("Code", "Function", "plotting_functions.R"))


#Note most plots are touched up in illustator because
#the legends are being difficult. However, the raw forms 
#should look 95% like the given figures 

#-------------------------------------------------------
# FIGURE 2 - Explore the effect of vector control on RE
#-------------------------------------------------------
# Simulate the model
source("Code/01_RE_mortality_P.R")

# Plot the composite 9-panel (Figure 2).",

# Retrieve "standard", "no_diff","better_m"
RE_mortality_P_SB <- subset(
  RE_mortality_P_post,
  RE_mortality_P_post$param %in% c("standard", "no_diff", "better_m")
)

# For better labeling of figure
RE_mortality_P_SB$param <- factor(RE_mortality_P_SB$param,
  levels = c("standard", "no_diff", "better_m"),
  labels = c("Standard", "Same", "Higher")
)

# Get the disturbance time
dstb_time <- get_parameters("standard")["disturbance_time"]

### Ensure that the RE legend is the same across different scenarios
RE_limits <- c(
  round(min(RE_mortality_P_SB$RE), 1),
  round(max(RE_mortality_P_SB$RE), 1) + 0.1
)

# Plotting RE over time
Panel_RE_Time <-
  plot_RE_dynamics(RE_mortality_P_SB, dstb_time, "No", NA) +
  facet_wrap(~param, ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11)
  )

# Plot total vector abundance over time with the RE As color
Panel_Vec_Time <- plot_NV_RE(RE_mortality_P_SB, dstb_time, "No", RE_limits) +
  facet_wrap(~param, ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

# Remove the situation where 100% of the primary vector is removed,
# not interesting (dynamically, and a little bug that occurs in Panel B!
removed_0 <- subset(RE_mortality_P_SB, RE_mortality_P_SB$id != 0)

# Plot the secondary versus primary vector with the RE as color
Panel_Vec_Time <- plot_NP_NM_RE(removed_0, postdisturb = "No", RE_limits) +
  facet_wrap(~param, ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

# Full composite figure with all panels together in a 3x3
Panel_RE_Time /
  (Panel_Vec_Time + theme(aspect.ratio = 1)) /
  (Panel_Vec_Time + theme(aspect.ratio = 1)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


ggsave(here("Figures_Process", "Figure_2.pdf"),
  width = 6, height = 6.5, units = "in"
)

#--------------------------------------------------
# FIGURE 3 - Plot the S. Vector contribution to RE
#--------------------------------------------------
# Simulate the model
source("Code/02_M_Vec_Contribution.R")

plot_m_contribution_heatmap(RE_SM)[[1]] +
  theme(legend.position = "top") +
  coord_equal() +
  plot_m_contribution_lineplot(RE_SM_2_subset, dstb_time)

ggsave(here("Figures_Process", "Figure_3.pdf"),
  width = 8, height = 6, units = "in"
)


#---------------------------------------------------------
# FIGURE 4 - Plot the effect of interspecific competition
#---------------------------------------------------------
source("Code/03_Interspecific_Competition.R")

# The main competition plot
competition_plot <- plot_competition_PM(RE_COMPETITION)[[1]]

# Underlying phaseplots
combined_phaseplot <- (plot_phaseplot(df_075, isocline_df[[1]]) +
  theme(legend.position = "bottom")) +
  (plot_phaseplot(df_1, isocline_df[[2]]) + theme(legend.position = "bottom")) +
  (plot_phaseplot(df_125, isocline_df[[3]]) + theme(legend.position = "bottom"))

competition_plot
ggsave(here("Figures_Process", "Figure_4_A.pdf"),
  width = 6, height = 3, units = "in"
)

combined_phaseplot
ggsave(here("Figures_Process", "Figure_4_B.pdf"),
  width = 6, height = 3, units = "in"
)


#----------------------------------------------------------------------------
# FIGURE 5 - Plot the difference in biting rate and transmission probability
#----------------------------------------------------------------------------

# Model and plot (because the plotting is simple) are all in the same file
source("Code/04_secondary_ecology.R")

#-----------------------------------------
### Supplementary Material
#-----------------------------------------
#----------------------------------------------
# Plot RE for different disturbance intensities
#----------------------------------------------
source("Code/Supplementary Code/SUPP_Intensity_PM.R")

plot_heatmapintensity(df_expand_RE)

ggsave(here("Figures_Process", "R0_heatmap.pdf"),
       width = 9, height = 7, units = "in"
)

#----------------------------------------------
# Plot to show that increase in RE is due to 
# secondary vector
#----------------------------------------------
source("Code/Supplementary Code/SUPP_Abundance_dynamics.R")

#-------------------------------------------------
# Plot the composite three-panel Figure 1 but for
# the other cases
#-------------------------------------------------
# Retrieve "standard" and "better_m"
RE_mortality_P_worse <- subset(
  RE_mortality_P_post,
  RE_mortality_P_post$param %in%
    c("worse_m", "nonesec")
)

RE_mortality_P_worse$param <- factor(RE_mortality_P_worse$param,
  levels = c("worse_m", "nonesec"),
  labels = c("Worse", "None")
)

### Makes it easyy
RE_limits_worse <- c(
  round(min(RE_mortality_P_worse$RE), 1),
  round(max(RE_mortality_P_worse$RE), 1) + 0.1
)

Panel_RE_Time_SUPP <-
  plot_RE_dynamics(RE_mortality_P_worse, dstb_time, "No", NA) +
  facet_wrap(~param, ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11)
  )

# Plot total vector abundance over time with the RE As color
Panel_A_SUPP <- plot_NV_RE(RE_mortality_P_worse, dstb_time, "No", RE_limits) +
  facet_wrap(~param, ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

# Remove the situation where 100% of the primary vector is removed,
# not interesting (dynamically, and a little bug that occurs in Panel B!
# If you do it individually, it works :/)
removed_0 <- subset(RE_mortality_P_worse, RE_mortality_P_worse$id != 0)

# Plot the secondary versus primary vector with the RE as color
Panel_B_SUPP <- plot_NP_NM_RE(removed_0, postdisturb = "No", RE_limits) +
  facet_wrap(~param, ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

# Full composite figure with all
Panel_RE_Time_SUPP /
  (Panel_A_SUPP + theme(aspect.ratio = 1)) /
  (Panel_B_SUPP + theme(aspect.ratio = 1)) +
  plot_layout(guides = "collect")


ggsave(here("Figures_Process", "SUPP_Figure_1.pdf"),
  width = 6, height = 7, units = "in"
)

#-------------------------------------------------
# Plot how RE changes with c_MM
#-------------------------------------------------
source("Code/Supplementary Code/SUPP_Deviation_c_MM.R")

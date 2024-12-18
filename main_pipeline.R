library(makepipe)
library(here)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Rcpp)
library(viridis)
library(patchwork)
source(here("Code", "Function", "simulate_model_Functions.R"))
source(here("Code", "Function", "calculations.R"))
sourceCpp(here("Code", "Function", "onepatch_model.cpp"))

make_with_source(
  note = "Simulate model with different disturbance intensity for primary vector",
  source = "Code/Simulate_MortP.R",
  targets = "Output/RE_mortality_P_post.rds")

make_with_recipe(
  note = "Plot the RE the different disturbance intensity for primary vector",
  label = "plot_NP_NM_RE",
  recipe = {
    RE_mortality_P_post <- readRDS("Output/RE_mortality_P_post.rds")
    plot_NP_NM_RE(RE_mortality_P_post)
    ggsave(here("Figures_Process", "NP_NM_RE.pdf"),
      width = 7, height = 6, units = "in"
    )
  },
  targets = "Figures_Process/NP_NM_RE.pdf",
  dependencies = "Output/RE_mortality_P_post.rds",
  envir = environment()
)

make_with_recipe(
  note = "Plot the dynamics of the states for when disturbance intensity is 25% ",
  label = "plot_state_dynamics",
  recipe = {
    RE_mortality_P_post <- readRDS("Output/RE_mortality_P_post.rds")
    plot_state_dynamics(RE_mortality_P_post)
    ggsave(here("Figures_Process", "states_dynamics_GG.pdf"),
      width = 4, height = 8, units = "in"
    )
  },
  targets = "Figures_Process/states_dynamics_GG.pdf",
  dependencies = "Output/RE_mortality_P_post.rds",
  envir = environment()
)

make_with_recipe(
  note = "Plot main figure (Figure 1) which is a composite",
  label = "plot_figure_1",
  recipe = {
    RE_mortality_P_post <- readRDS("Output/RE_mortality_P_post.rds")
    plot_NP_NM_RE(RE_mortality_P_post) +
      plot_state_dynamics(RE_mortality_P_post)

    ggsave(here("Figures_Process", "Figure1.pdf"),
      width = 14, height = 7, units = "in"
    )
  },
  targets = "Figures_Process/Figure1.pdf",
  dependencies = "Output/RE_mortality_P_post.rds",
  envir = environment()
)

make_with_source(
  note = "Calculate RE (or R0) for different vector abundance",
  source = "Code/Simulate_R0.R",
  targets = "Output/df_expand_RE.rds"
)

make_with_recipe(
  note = "Plot R0 for different vector abundance",
  label = "plot_heatmap_R0",
  recipe = {
    df_expand_RE <- readRDS("Output/df_expand_RE.rds")
    plot_heatmapR0(df_expand_RE )
    ggsave(here("Figures_Process", "R0_heatmap.pdf"),
      width = 9, height = 7, units = "in"
    )
  },
  targets = "Figures_Process/R0_heatmap.pdf",
  dependencies = "Output/df_expand_RE.rds",
  envir = environment()
)

make_with_source(
  note = "Calculate RE (or R0) for different disturbance intensity",
  source = "Code/Simulate_Intensity_PM.R",
  targets = "Output/Maximum_RE_Mort_PM.rds"
)

make_with_recipe(
  note = "Plot RE for different disturbance intensities",
  label = "plot_heatmap_disturbance_intensities",
  recipe = {
    df_expand_RE <- readRDS("Output/Maximum_RE_Mort_PM.rds")
    plot_heatmapR0(df_expand_RE)
    ggsave(here("Figures_Process", "R0_heatmap.pdf"),
      width = 9, height = 7, units = "in"
    )
  },
  targets = "Figures_Process/Maximum_RE_Mort_PM.pdf",
  dependencies = "Output/df_expand_RE.rds",
  envir = environment()
)


###Supplementary Material
make_with_recipe(
  note = "Plot supplementary material to compare vector abundance
  and maximum RE",
  label = "plot_supp_1",
  recipe = {
    RE_mortality_P_post <- readRDS("Output/RE_mortality_P_post.rds")
    plot_comparison_RE(RE_mortality_P_post, "id")

    ggsave(here("Figures_Process", "Supp1.pdf"),
      width = 9, height = 7, units = "in")
  },
  targets = "Figures_Process/Supp1.pdf",
  dependencies = "Output/RE_mortality_P_post.rds",
  envir = environment()
)


show_pipeline()

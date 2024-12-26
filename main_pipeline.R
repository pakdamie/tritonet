source(here("Code", "Function", "packages.R"))
sourceCpp(here("Code", "Function", "onepatch_model.cpp"))
source(here("Code", "Function", "simulate_functions.R"))
source(here("Code", "Function", "calculate_functions.R"))
source(here("Code", "Function", "plot_functions.R"))

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
      width = 7, height = 6, units = "in")
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
    plot_heatmapR0(df_expand_RE)
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
  note = "Plot the vector abundance for the maximum RE and the maximum
  vector abundance and corresponding RE",
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

make_with_recipe(
  note = "Plot total vector abundance over time and highlight when 
  the maximum RE is reached",
  label = "plot_totalV_maxRE",
  recipe = {
    RE_mortality_P_post <- readRDS("Output/RE_mortality_P_post.rds")
    plot_NV_RE(RE_mortality_P_post)
    
    ggsave(here("Figures_Process", "plot_totalV_maxRE.pdf"),
           width = 9, height = 7, units = "in")
  },
  targets = "Figures_Process/plot_totalV_maxRE.pdf",
  dependencies = "Output/RE_mortality_P_post.rds",
  envir = environment()
)

make_with_recipe(
  note = "Plot total vector abundance over time and highlight when 
  the maximum RE is reached (assuming vectors are the same)",
  label = "plot_totalV_maxRE_same",
  recipe = {
    RE_mortality_P_same <- readRDS("Output/RE_mortality_P_same_DF.rds")
    plot_state_dynamics(RE_mortality_P_same, mortality_P = 0.01)
    
    
    plot_NV_RE(RE_mortality_P_same) + scale_color_viridis()
    plot_comparison_RE(RE_mortality_P_same, "id")
    ggsave(here("Figures_Process", "plot_totalV_maxRE_same.pdf"),
           width = 9, height = 7, units = "in")
  },
  targets = "Figures_Process/plot_totalV_maxRE_same.pdf",
  dependencies = "Output/RE_mortality_P_same_DF.rds",
  envir = environment()
)



show_pipeline()

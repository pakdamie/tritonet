#This is the pipeline for the entirety of the manuscript.
# Maintained by Damie Pak
library(here)
source(here("Code", "Function", "packages.R"))
sourceCpp(here("Code", "Function", "model_vectors_host.cpp"))
source(here("Code", "Function", "simulate_functions.R"))
source(here("Code", "Function", "calculate_functions.R"))
source(here("Code", "Function", "plotting_functions.R"))

make_with_source(
  note = "Simulate model with different disturbance intensity for p. vector",
  source = "Code/Simulate_MortP.R",
  targets = "Output/RE_mortality_P_post.rds")

make_with_recipe(
  note = "Plot the composite three-panel Figure 1.",
  label = "plot_NP_NM_RE",
  recipe = {
    RE_mortality_P_post <- readRDS("Output/df_RE_mortality_P_R0.rds")

    ### Makes it easyy
    RE_limits <- c(
      round(min(RE_mortality_P_post$RE), 1),
      round(max(RE_mortality_P_post$RE), 1) + 0.1
    )

    Panel_RE_Time <- plot_RE_dynamics(RE_mortality_P_post, "No", NA)

    # Plot total vector abundance over time with the RE As color
    Panel_A <- plot_NV_RE(RE_mortality_P_post, "No", RE_limits)

    # Remove the situation where 100% of the primary vector is removed,
    # not interesting (dynamically, and a little bug that occurs in Panel B!
    # If you do it individually, it works :/)
    removed_0 <- subset(RE_mortality_P_post, RE_mortality_P_post$id != 0)

    # Plot the secondary versus primary vector with the RE as color
    Panel_B <- plot_NP_NM_RE(removed_0, postdisturb = "No", RE_limits)

    layout <- c(
      area(1, 1, 2, 2),
      area(1, 3, 1),
      area(2, 3, 2, 3)
    )


    # Full composite figure with all
    Panel_RE_Time + Panel_A + theme(aspect.ratio = 1) + 
      Panel_B + theme(aspect.ratio = 1) +
      plot_layout(guides = "collect", design = layout) +
      plot_annotation(
        tag_levels = c("A", "1"),
        tag_sep = ".", tag_suffix = "."
      )

    ggsave(here("Figures_Process", "Figure_1.pdf"),
      width = 10, height = 8, units = "in"
    )
  },
  targets = "Figures_Process/Figure_1.pdf",
  dependencies = "Output/RE_mortality_P_post.rds",
  envir = environment()
)


make_with_recipe(
  note = "Plot the composite three-panel Figure 1.",
  label = "plot_NP_NM_RE",
  recipe = {
# Remove the situation where 100% of the primary vector is removed,
# not interesting (dynamically, and a little bug that occurs in Panel B!
# If you do it individually, it works :/)
removed_0 <- subset(RE_mortality_P_post, RE_mortality_P_post$id != 0)

plot_NM_REff(removed_0,"No") + 

}


make_with_source(
  note = "Calculate RE (or R0) for different vector abundance",
  source = "Code/Simulate_R0.R",
  targets = "Output/df_expand_RE.rds"
)

make_with_recipe(
  note = "Plot R0 for different vector abundance using different
  parameters",
  label = "plot_heatmap_R0",
  recipe = {
    df_expand_RE <- readRDS("Output/df_expand_RE.rds")
    plot_heatmapR0(df_expand_RE)
    ggsave(here("Figures", "R0_heatmap.pdf"),
      width = 13, height = 5, units = "in"
    )
    
  },
  targets = "Figures/R0_heatmap.pdf",
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
      width = 7, height = 6, units = "in")
  },
  targets = "Figures_Process/Supp1.pdf",
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
    
    
    plot_NV_RE(RE_mortality_P_same, "Yes") + scale_color_viridis()
    plot_comparison_RE(RE_mortality_P_same, "id")
    ggsave(here("Figures_Process", "plot_totalV_maxRE_same.pdf"),
           width = 9, height = 7, units = "in")
  },
  targets = "Figures_Process/plot_totalV_maxRE_same.pdf",
  dependencies = "Output/RE_mortality_P_same_DF.rds",
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


show_pipeline()

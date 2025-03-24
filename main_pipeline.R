#This is the pipeline for the entirety of the manuscript.
# Maintained by Damie Pak

library(here)
source(here("Code", "Function", "packages.R"))
sourceCpp(here("Code", "Function", "model_vectors_host.cpp"))
source(here("Code", "Function", "simulate_functions.R"))
source(here("Code", "Function", "calculate_functions.R"))
source(here("Code", "Function", "plotting_functions.R"))

#---------
#FIGURE 2
#---------
make_with_source(
  note = "Simulate model with different disturbance intensity for p. vector 
  and different parameter",
  source = "Code/01_RE_mortality_P.R",
  targets = "Output/RE_mortality_P_post.rds")

make_with_recipe(
  note = "Plot the composite 9-panel (Figure 2).",
  label = "plot_NP_NM_RE",
  recipe = {
    
    RE_mortality_P_post <- readRDS("Output/df_RE_mortality_P_R0.rds")

    #Retrieve "standard", "no_diff","better_m"
    
    RE_mortality_P_SB <- subset(RE_mortality_P_post, 
                                RE_mortality_P_post$param %in%
                                c("standard", "no_diff","better_m"))
    
    RE_mortality_P_SB$param <- factor(RE_mortality_P_SB$param,
                                         levels = c("standard", "no_diff", "better_m"),
                                         labels = c("Standard", "Same", "Higher"))
    
    dstb_time<- get_parameters("standard")["disturbance_time"]
    
    ### Makes it easyy
    RE_limits <- c(
      round(min(RE_mortality_P_SB$RE), 1),
      round(max(RE_mortality_P_SB$RE), 1) + 0.1
    )

    Panel_RE_Time <- 
      plot_RE_dynamics(RE_mortality_P_SB,  dstb_time, "No", NA) + 
      facet_wrap(~param, ncol = 3) +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = 11))

    # Plot total vector abundance over time with the RE As color
    Panel_A <- plot_NV_RE(RE_mortality_P_SB,  dstb_time,"No", RE_limits) + 
      facet_wrap(~param,ncol = 3) +
      theme(strip.background = element_blank(),
            strip.text = element_blank()) 

    # Remove the situation where 100% of the primary vector is removed,
    # not interesting (dynamically, and a little bug that occurs in Panel B!
    # If you do it individually, it works :/)
    removed_0 <- subset(RE_mortality_P_SB, RE_mortality_P_SB$id != 0)

    # Plot the secondary versus primary vector with the RE as color
    Panel_B <- plot_NP_NM_RE(removed_0, postdisturb = "No", RE_limits)+ 
      facet_wrap(~param,ncol = 3) +
      theme(strip.background = element_blank(),
            strip.text = element_blank())

  
    # Full composite figure with all
    Panel_RE_Time  / 
      (Panel_A + theme(aspect.ratio = 1)) / 
      (Panel_B + theme(aspect.ratio = 1)) + 
      plot_layout(guides = "collect") &
      theme(legend.position = 'bottom')
    
    
    ggsave(here("Figures_Process", "Figure_2.pdf"),
      width = 6, height =6.5, units = "in"
    )
  },
  targets = "Figures_Process/Figure_2.pdf",
  dependencies = "Output/df_RE_mortality_P_R0.rds",
  "Code/Functions/plotting_functions.R",
  envir = environment()
)

#---------
#FIGURE 3
#---------
make_with_source(
  note = "Simulate model to investigate secondary contribution",
  source = "Code/02_M_Vec_Contribution.R",
  targets = "Output/RE_SM_2_subset.rds", "Output/RE_SM.rds")


make_with_recipe(
  note = "Plot the S. Vector contribution to RE (Figure 3)",
  label = "plot_M_RE",
  recipe = {
    
    dstb_time<- get_parameters("standard")["disturbance_time"]
    
    
    RE_SM <- readRDS("Output/RE_SM.rds") #Non-temporal
    RE_SM_2_subset <- readRDS("Output/RE_SM_2_subset.rds") #Temporal
    
    (plot_m_contribution_heatmap (RE_SM)[[1]] + theme(legend.position = 'top')+
        coord_equal()) +
      plot_m_contribution_lineplot(RE_SM_2_subset,dstb_time)
    
    ggsave(here("Figures_Process", "Figure_3.pdf"),
           width = 8, height = 6, units = "in"
    )
  },
  targets = "Figures_Process/Figure_3.pdf",
  dependencies = "Output/RE_SM.rds","Output/RE_SM_2_subset.rds",
  envir = environment()
)

#---------
#FIGURE 4
#---------


# Remove legend from competition plot
competition_plot <- plot_competition_PM(RE_COMPETITION)[[1]] 
# Retain legends for phaseplots
combined_phaseplot <- (plot_phaseplot(df_075, isocline_df[[1]]) + 
                         theme(legend.position = 'bottom')) +
  (plot_phaseplot(df_1, isocline_df[[2]]) + theme(legend.position = 'bottom')) +
  (plot_phaseplot(df_125, isocline_df[[3]]) + theme(legend.position = 'bottom')) 




competition_plot
 ggsave(here("Figures_Process", "Figure_4_A.pdf"),
        width = 6, height = 3, units = "in"
 )
 combined_phaseplot 
ggsave(here("Figures_Process","Figure_4_B.pdf"),
                    width = 6, height = 3, units = "in")


plot_competition_PM(RE_COMPETITION)[[2]] 
ggsave(here("Figures_Process", "Supp_Secondary_Vectors.pdf"),
       width = 6, height = 6, units = "in"
)















###Supplementary Material

make_with_source(
  note = "Plot the abundance of NP/NM over time",
  source = "Code/Supplementary Code/SUPP_Abundance_dynamics.R",
  targets = "Figures/SUPP_abundance_dynamics.pdf"
)



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

make_with_recipe(
  note = "Plot the composite three-panel Figure 1 but for
  the other standards",
  label = "SUPP_plot_NP_NM_RE",
  recipe = {
    
    RE_mortality_P_post <- readRDS("Output/df_RE_mortality_P_R0.rds")
    
    dstb_time<- get_parameters("standard")["disturbance_time"]
    
    
    #Retrieve "standard" and "better_m"
    RE_mortality_P_worse <- subset(RE_mortality_P_post, 
                                   RE_mortality_P_post$param %in%
                                   c("worse_m", "nonesec"))
    
    RE_mortality_P_worse$param <- factor(RE_mortality_P_worse$param,
                                      levels = c("worse_m","nonesec"),
                                      labels = c("Worse", "None"))
    
    ### Makes it easyy
    RE_limits_worse <- c(
      round(min(RE_mortality_P_worse$RE), 1),
      round(max(RE_mortality_P_worse$RE), 1) + 0.1
    )
    
    Panel_RE_Time_SUPP <- 
      plot_RE_dynamics(RE_mortality_P_worse ,dstb_time, "No", NA) + 
      facet_wrap(~param, ncol = 3) +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = 11))
    
    # Plot total vector abundance over time with the RE As color
    Panel_A_SUPP <- plot_NV_RE(RE_mortality_P_worse,dstb_time, "No", RE_limits) + 
      facet_wrap(~param,ncol = 3) +
      theme(strip.background = element_blank(),
            strip.text = element_blank()) 
    
    # Remove the situation where 100% of the primary vector is removed,
    # not interesting (dynamically, and a little bug that occurs in Panel B!
    # If you do it individually, it works :/)
    removed_0 <- subset(RE_mortality_P_worse,RE_mortality_P_worse$id != 0)
    
    # Plot the secondary versus primary vector with the RE as color
    Panel_B_SUPP <- plot_NP_NM_RE(removed_0, postdisturb = "No", RE_limits)+ 
      facet_wrap(~param,ncol = 3) +
      theme(strip.background = element_blank(),
            strip.text = element_blank())
    
    
    # Full composite figure with all
    Panel_RE_Time_SUPP  / 
      (Panel_A_SUPP + theme(aspect.ratio = 1)) / 
      (Panel_B_SUPP + theme(aspect.ratio = 1)) + 
      plot_layout(guides = "collect") 
    
    
    ggsave(here("Figures_Process", "SUPP_Figure_1.pdf"),
           width = 6, height = 7, units = "in"
    )
  },
  targets = "Figures_Process/SUPP_Figure_1.pdf",
  dependencies = "Output/RE_mortality_P_post.rds",
  envir = environment()
)










show_pipeline()

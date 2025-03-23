# --------------------------------------------------------------------------------
# Research Question:
# How does the change in the primary to secondary coefficient influence RM
# with the change in the transmission efficiency?
# --------------------------------------------------------------------------------
param_standard <- get_parameters("standard")
c_PM_standard <- param_standard["c_PM"] ## Competition effect of p.vector on s.vector
c_MM_standard <- param_standard["c_MM"]
c_MP_standard <- param_standard["c_MP"]
c_PP_standard <- param_standard["c_PP"]

modifier <- c(0.75,1,1.25,seq(.6766667,2,0.01))
Mortality_P <- c(0.01,0.25, 0.5,0.75)


competition_param <-
  data.frame(expand.grid(
    c_PM = c_PM_standard * modifier,
    mortality_P = Mortality_P
  ))


competition_param_list <- vary_parameter_value(
  param_standard, c( "c_PM", "mortality_P"), competition_param
)

# --------------------------------------------------------------------------------
# Step 1: Simulate Model Output and Calculate Re
# --------------------------------------------------------------------------------



RE_COMPETITION <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("c_PM", "mortality_P"),
    vector_value = competition_param
  ) |>
  Calculate_change_baseline(
    competition_param_list,
    competition_param, "No"
  )

#Note: When you analytically solve for when the primary vector competitively 
#excludes the secondary vector it should be when C_PM is multiplied by 1.5.
#There is a slight numerical error when simulating it so I decided to fix it myself
#It's less obvious when you simulate for a loooooong time, but it's a lot of
#computational load. 

RE_COMPETITION$standardized <- RE_COMPETITION$c_PM/c_PM_standard
RE_COMPETITION[RE_COMPETITION$standardized > 1.49,]$RE <- 0
RE_COMPETITION[RE_COMPETITION$standardized > 1.49,]$max_NM <- 0

# --------------------------------------------------------------------------------
# Step 2: Save Processed Data
# --------------------------------------------------------------------------------

# Save the filtered dataset as an RDS file in the "Output" directory

saveRDS(RE_COMPETITION, file = here("Output", "RE_COMPETITION.rds"))

# --------------------------------------------------------------------------------
# Research Question:
# How does the change in the primary to secondary coefficient influence RM
# with the change in the transmission efficiency?
# --------------------------------------------------------------------------------


scaling_factors <- c(0.75,1,1.25)

isocline_df <- lapply(scaling_factors, function(factor) {
  list(
    Isocline = Calculate_Isocline(factor * c_PM_standard),
    Phaseplot = Calculate_Phaseplot(factor * c_PM_standard)
  )
})

competition_param_001<-
  data.frame(expand.grid(
    c_PM = c_PM_standard * scaling_factors ,
    mortality_P = 0.01
  ))


competition_param_list_001 <- vary_parameter_value(
  param_standard, c( "c_PM", "mortality_P"), 
  competition_param_001
)

RE_COMPETITION_2 <-
  Simulate_Model_Output(
    parameter = get_parameters("standard"),
    infection_start = "No",
    variable_interest = c("c_PM", "mortality_P"),
    vector_value = competition_param_001
  ) 



process_testing <- function(RE_data, params) {
  testing <- Calculate_Human_Reff_Expanded(RE_data, params)
  df_interest<- testing[testing$time > params[["disturbance_time"]] - 1, ]
  return(df_interest)
}

df_075 <- Process_Results(RE_COMPETITION_2[[1]],competition_param_list[[1]] )
df_1 <- Process_Results(RE_COMPETITION_2[[2]],competition_param_list[[2]] )
df_125 <- Process_Results(RE_COMPETITION_2[[3]],competition_param_list[[3]] )



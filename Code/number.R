### Run this script for numerical results for 
min_RE_fig_1  <- min(subset(RE_mortality_P_post, RE_mortality_P_post$id == 0.01)$RE); 
paste("The minimum RE reached in Figure 1 is:", round(min_RE_fig_1,2), sep = " ")

max_RE_fig_1  <- max(subset(RE_mortality_P_post, RE_mortality_P_post$id == 0.01)$RE); 
paste("The maximum RE reached in Figure 1 is:", round(max_RE_fig_1,2), sep = " ")

max_NM_00_fig_1  <- max(subset(RE_mortality_P_post, RE_mortality_P_post$id %in% 0)$NM); 
max_NM_01_fig_1  <- max(subset(RE_mortality_P_post, RE_mortality_P_post$id %in% 0.01)$NM); 
max_NM_25_fig_1  <- max(subset(RE_mortality_P_post, RE_mortality_P_post$id %in% 0.25)$NM); 
eq_NM_25_fig_1  <- (subset(RE_mortality_P_post, RE_mortality_P_post$time == 9124)$NM); 


paste("The maximum NM reached in Figure 1 is:", round(max_NM_01_fig_1,2), sep = " ")
paste("The maximum NM reached in Figure 1 is:", round(max_NM_25_fig_1,2), sep = " ")

max(RE_SECONDARY$RE)

calculate_predisturb_initial_list <- function(){
        
  param_no_disturb <- get_parameters("no_disturb")
  Initial_List<- create_initial_states(get_parameters("standard"))
  Mod_Predisturb <- discrete_trito_model_rcpp_ONEPATCH(
        HS = Initial_List[[1]],
        HI = Initial_List[[2]],
        HR = Initial_List[[3]],
        PS = Initial_List[[4]],
        PI = Initial_List[[5]],
        MS = Initial_List[[6]],
        MI = Initial_List[[7]],
        param = param_no_disturb)

  Predisturb_EqStates<- lapply(Mod_Predisturb,
                             function(x) x[((365 * 25)),])

  eq_initial_list<- create_initial_states(get_parameters("post_disturb"))
  eq_initial_list[[1]][1] <- Predisturb_EqStates[[1]]
  eq_initial_list[[2]][1] <- Predisturb_EqStates[[2]]
  eq_initial_list[[3]][1] <- Predisturb_EqStates[[3]]
  eq_initial_list[[4]][1] <- Predisturb_EqStates[[4]]
  eq_initial_list[[5]][1] <- Predisturb_EqStates[[5]]
  eq_initial_list[[6]][1] <- Predisturb_EqStates[[6]]
  eq_initial_list[[7]][1] <- Predisturb_EqStates[[7]]
  
  return(eq_initial_list)
}
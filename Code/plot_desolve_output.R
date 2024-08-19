plotter_dynamics <- function(desolve_list) {
        
  results <- desolve_list[[1]]
  chosen_patch <- desolve_list[[2]]
  
  ### TOTAL HUMAN- N_H
  human_all <- cbind(
    time = results[, "time"],
    results[, grepl("HS", names(results)) |
      grepl("HI", names(results)) |
      grepl("HR", names(results))]
  )

  ### Total Human-Infected
  human_infected <- cbind(
    time = results[, "time"],
    results[, grepl("HI", names(results))]
  )

  ### Total susceptible vectors both primary and secondary
  vector_sus <- cbind(
    time = results[, "time"],
    results[grepl("PS", names(results)) |
      grepl("SS", names(results))]
  )
  ### Total infected vectors both primary and secondary
  vector_inf <- cbind(
    time = results[, "time"],
    results[grepl("PI", names(results)) |
      grepl("SI", names(results))]
  )
  ### Total primary vectors (susceptible and infected)
  vector_P <- cbind(
    time = results[, "time"],
    results[grepl("PS", names(results)) |
      grepl("PI", names(results))]
  )

  ### Total secondary vectors (susceptible and infected)
  vector_S <- cbind(
    time = results[, "time"],
    results[grepl("SS", names(results)) |
      grepl("SI", names(results))]
  )

  vector_P_sus_only <- cbind(
          time = results[, "time"],
          results[grepl("PS", names(results))])
  
  vector_P_inf_only <- cbind(
          time = results[, "time"],
          results[grepl("PI", names(results))])
  
  vector_P_all_patch <- cbind(time = results[, "time"],
                              vector_P_sus_only[,2:ncol( vector_P_sus_only)]+
                         vector_P_inf_only[,2:ncol( vector_P_inf_only)])
  colnames( vector_P_all_patch ) <- c("time", paste0("P", seq(1,ncol( vector_P_inf_only)-1)))

  
  vector_S_sus_only <- cbind(
          time = results[, "time"],
          results[grepl("SS", names(results))])
  
  vector_S_inf_only <- cbind(
          time = results[, "time"],
          results[grepl("SI", names(results))])
  
  vector_S_all_patch <- cbind(time = results[, "time"],
                              vector_S_sus_only[,2:ncol( vector_S_sus_only)]+
                                      vector_S_inf_only[,2:ncol( vector_S_inf_only)])
  colnames( vector_S_all_patch ) <- c("time", paste0("S", seq(1,ncol( vector_S_inf_only)-1)))
  
  
  ### A function just to melt and label which patch
  ### has been sprayed
  melt_label_subpop <- function(df, group_name, chosen_patch) {
    melted_df <- melt(df, id.vars = "time")
    melted_df$sprayed <- "NO"

    if (group_name %in% c("HR", "HI", "HR")) {
      patch_number <- parse_number(chosen_patch)
      patch_sprayed <- paste0(group_name, patch_number)
    } else {
      patch_sprayed <- chosen_patch
    }

    melted_df$sprayed[melted_df$variable %in% patch_sprayed] <- "YES"
    return(melted_df)
  }
  plot_lineplot <- function(melted_df, title, facet_group) {
 melted_df$grouping <-  sub("^([[:alpha:]]*).*", "\\1", melted_df$variable)
 
 melted_df$number <- parse_number(as.character(melted_df$variable))
    base_GG <- ggplot(melted_df, aes(x = time, y = log10(value+1), group = variable)) +
      xlab("Time") +
      ylab("Abundance") +
      ggtitle(title) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15)
      )
    
    if(facet_group == "yes"){

            base_GG2 =  base_GG+ 
                    geom_line(aes(linetype = sprayed,color = number ))+ 
                    scale_color_viridis(option = "viridis", discrete = FALSE) +
                    facet_wrap(~grouping)
    }
    else{
            base_GG2 = base_GG + geom_line(aes(linetype = sprayed,
                                               color = number ))+
                    scale_color_viridis(option = "viridis", discrete = FALSE) 
                    
    }
          return(base_GG2)
  }
  
  
  melted_H_sus <- melt_label_subpop(  human_all, "HS",
                                      chosen_patch)
  
  melted_H_infected <- melt_label_subpop(human_infected, "HI", chosen_patch)
  melted_vec_infected <- melt_label_subpop(  vector_inf, NA, chosen_patch)
  melted_vec_susceptible <- melt_label_subpop(  vector_sus, NA, chosen_patch)
  
  
  melted_p_all <-  melt_label_subpop(  vector_P_all_patch , NA, chosen_patch)
  melted_s_all <-  melt_label_subpop(  vector_S_all_patch , NA, chosen_patch)
  
  GG_humans_all <- plot_lineplot( melted_H_sus , "All Human","yes") + 
          facet_wrap(~fct_relevel(grouping,'HS','HI','HR'))
  GG_vectors_inf <- plot_lineplot( melted_vec_infected , "Infected vectors","yes")
  GG_vectors_sus <- plot_lineplot( melted_vec_susceptible, "Sus. vectors","yes")
  
  GG_vectors_P <- plot_lineplot(   melted_p_all  , "All primary","no")
  GG_vectors_S <- plot_lineplot(   melted_s_all , "All secondary","no")
  
  return(list(GG_humans_all,
              GG_vectors_inf,
              GG_vectors_sus,
              GG_vectors_P + GG_vectors_S ))
  
}

(plotter_dynamics(results_ode)[[4]][[1]]+ylim(-1,3))+
      (plotter_dynamics(results_ode)[[4]][[2]]+ylim(-1,3))
        


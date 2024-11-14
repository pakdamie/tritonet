#This script is for simulating a lot of infections
#based on different connectance, spatial coverage, 
#and impact on primary.
impact_on_primary <- seq(0,1,0.10)
spatial_coverage <- seq(0,1,0.10)

expanded_param <- expand.grid(
  impact_P = impact_on_primary, 
  sp_cov = spatial_coverage,
  connectance = seq(1,10),
  replicate = seq(1,100))


###THE PARAMETERS THAT I'M INTERESTED IN

###The function to simulate a lot of infection
simulate_infection <- function(params) {
  impact_P <- params$impact_P
  sp_cov <- params$sp_cov
  connectance <- params$connectance
  replicate <- params$replicate
  
  # Call your discrete_trito_model with the correct parameters
  return(
    discrete_trito_model(1000, param, 500, impact_P, 0.9, sp_cov, 
                              igraph_list[[connectance]][[replicate]], 1)
    )
}








#Batch processing
for (batch in 1:num_batches) {
  # Determine start and end indices for the current batch
  start <- (batch - 1) * batch_size + 1
  end <- min(batch * batch_size, nrow(expanded_param))
  
  # Run simulations in parallel for the current batch
  results_batch <- mclapply(start:end, function(i) {
    simulate_infection(expanded_param[i, ])
  }, mc.cores = num_cores - 2)
  
  # Save the results of the current batch
  save(results_batch, 
       file = here("Output", "Connectivity_Intensity", 
                   paste0("sim_result_", batch, ".RData")))
  
  # Remove the current batch results from memory
  rm(results_batch)
}
}

process_parallel 

num_files <- length(list.files(here("Output", "Connectivity_Intensity")))
data_list <- list()


for (i in 1:num_files) {
  # Dynamically create the file name
  file_name <- paste0("sim_result", i, ".RData")
  
  # Read in the file and store it in the list
   load(here("Output", "Connectivity_Intensity", file_name))
   
   calculated_RE <- NULL
   
   for (j in 1:100) {
    result = results_batch[[j]]
    calculated_RE [[j]] <-  calculate_R_effective_discrete_patch(param, result,500)
   }

data_list[[i]] <- do.call(rbind,calculated_RE)
}
data_list_2<- do.call(rbind,unlist(data_list, recursive = FALSE))


full_data_list <- cbind(data_list_2,expanded_param  )

###LET's look at 50% conenctivity

full_data_list50 <- subset(full_data_list, round(full_data_list$sp_cov,2) ==0.50)

splitted_full_data_list50 <- split(full_data_list50,list(full_data_list50$impact_P,full_data_list50$connectance))

appli<- lapply(splitted_full_data_list50,  function(x)
     cbind(CV = unique(max(x$CV)),RE=  unique(max(x$RE)), 
           impact_P = unique(x$impact_P),connectance= unique(x$connectance)))

connect_stuff <- do.call(rbind.data.frame, appli)


MAX_RE_50_GG <- ggplot(connect_stuff, 
  aes(x = connectance, y= 1-impact_P, fill = (RE))) + 
  geom_tile(color = 'black') +
  scale_x_continuous(expand = c(0,0), breaks = seq(1,10),
                     labels = seq(0.05, 0.5, 0.05))+
  scale_y_continuous(expand = c(0,0)) + 
  xlab("Connectance") + 
  ylab("Proportion of P. species removed ") + 
  scale_fill_viridis(option = 'viridis', name = "Maximum RE") + 
  ggtitle("90% spatial coverage") + 
  theme_classic()+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14.5))

  CV_RE_50_GG <- ggplot(connect_stuff, 
                         aes(x = connectance, y= 1-impact_P, fill = (CV))) + 
    geom_tile(color = 'black') +
    scale_x_continuous(expand = c(0,0), breaks = seq(1,10),
                       labels = seq(0.05, 0.5, 0.05))+
    scale_y_continuous(expand = c(0,0)) + 
    xlab("Connectance") + 
    ylab("Proportion of P. species removed ") + 
    scale_fill_viridis(option = 'viridis', name = "Maximum CV RE") + 
    ggtitle("90% spatial coverage") + 
  theme_classic() + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14.5))
  
  
MAX_RE_50_GG  +   CV_RE_50_GG 
  
ggsave(here("Figures_Process","max_cv_90_GG.pdf"), width = 13, height = 4, units = 'in')


###LET's look at 100% species removal

full_data_list100 <- subset(full_data_list, round(full_data_list$impact_P,2) ==0.70)

splitted_full_data_list100 <- split(full_data_list100,
                                    list(full_data_list100 $sp_cov,
                                         full_data_list100 $connectance))

appli100<- lapply(splitted_full_data_list100,  function(x)
  cbind(CV = unique(max(x$CV)),RE=  unique(max(x$RE)), 
        sp_cov = unique(x$sp_cov),connectance= unique(x$connectance)))

connect_stuff100 <- do.call(rbind.data.frame, appli100)


MAX_RE_100_GG <- ggplot(connect_stuff100, 
                       aes(x = connectance, y= sp_cov,
                           fill = (RE))) + 
  geom_tile(color = 'black') +
  scale_x_continuous(expand = c(0,0), breaks = seq(1,10),
                     labels = seq(0.05, 0.5, 0.05))+
  scale_y_continuous(expand = c(0,0)) + 
  xlab("Connectance") + 
  ylab("Spatial coverage") + 
  scale_fill_viridis(option = 'viridis', name = "Maximum RE") + 
  ggtitle("70% of the species removed") + 
  theme_classic()+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14.5))

CV_RE_100_GG <- ggplot(connect_stuff100, 
                      aes(x = connectance, y= sp_cov, fill = (CV))) + 
  geom_tile(color = 'black') +
  scale_x_continuous(expand = c(0,0), breaks = seq(1,10),
                     labels = seq(0.05, 0.5, 0.05))+
  scale_y_continuous(expand = c(0,0)) + 
  xlab("Connectance") + 
  ylab("Spatial coverage") + 
  scale_fill_viridis(option = 'viridis', name = "Maximum CV RE") + 
  ggtitle("70% of the species removed") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14.5))


MAX_RE_100_GG  +   CV_RE_100_GG 

ggsave(here("Figures_Process","max_cv_70_sp_GG.pdf"), width = 13, height = 4, units = 'in')

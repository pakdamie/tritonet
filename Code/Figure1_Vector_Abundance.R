
# Standard Parameter ---------------------------------------------------------------
param_standard <- c(
  b_H = 1/ (1000), ## Human mortality rate
  b_P = 0.01, # P. Vector birth rate
  b_M = 0.01, # S. Vector birth rate
  mu_H = 1/ (1000), ## Human death rate
  f_P = 0.02, # Biting rate of the p. vector
  f_M = 0.020 *0.75, # Biting rate of the s.vector
  theta_P = 0.70, # Transmission probability of p. vector
  theta_M = 0.70 * 0.75, # Transmission probability of s. vector
  theta_H = 0.50, # Transmission probability of human
  gamma = 1 / 90, # Recovery rate of infected human
  c_PM = 4e-6, ## Competition effect of p.vector on s.vector
  c_MP = 2e-6, ## Competition effect of s.vector on p.vector
  c_PP = 4.5e-6, ## Competition effect of p.vector on s.vector
  c_MM = 3e-6, ## Competition effect of s.vector on s.vector
  ntime = 365 * 50,
  disturbance_time = 365 * 25,
  delta_T = 1,
  prop = 1,
  mortality_P = 0.25, # This will change
  mortality_M = 1)







# Running models and cleaning it up ---------------------------------------
Mortality_P <- seq(0.05,1,0.05)
Mortality_P[1] <- 0.01

Model_mortality_P <- Simulate_Model_Output(
  param_standard, "mortality_P", Mortality_P)

RE_mortality_P <- lapply(Model_mortality_P, 
                         Calculate_analytical_REff,
                         param = param_standard)


#Assign mortality_P value

for (i in seq_along(RE_mortality_P)){
  RE_mortality_P[[i]]$id <- Mortality_P[[i]]        
}

RE_mortality_P_DF <- do.call(rbind, RE_mortality_P)
RE_mortality_P_post <- subset(RE_mortality_P_DF, 
                              RE_mortality_P_DF$time > 9120 &
                              RE_mortality_P_DF$time < 14000)




# Plotting ----------------------------------------------------------------


GG_PtoMabundance_RE <- ggplot(RE_mortality_P_post, 
       aes(x = NP, y = NM, color = RE)) +
  geom_path(linewidth = 1.2) +
  geom_point(data = subset(RE_mortality_P_post, 
                           RE_mortality_P_post$time == 9125),
             aes(x = NP, y= NM), color= 'black') +
  scale_colour_gradient2(
    low = "#448C81",
    mid = "#EFEAD7",
    high = "#C25D31",
    midpoint = 1,
    name = expression(R[E])) + 
  xlab(expression("Abundance of primary vectors " *"(" * N[P] * ")")) + 
  ylab(expression("Abundance of secondary vectors " *"(" * N[M] * ")")) +   
  theme_classic() + 
  theme(axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 15, color = 'black'))


Maximum_RE_DF <- 
  do.call(rbind,by(RE_mortality_P_post, 
  RE_mortality_P_post$id, 
  function(x) 
  x[which.max(x$RE),],
  simplify = TRUE))

Maximum_Abundance_DF <- 
  do.call(rbind, by(RE_mortality_P_post, 
  RE_mortality_P_post$id, 
  function(x) 
  x[which.max(x$NP + x$NM),],
  simplify = TRUE))

GG_Totalvectorabundance_RE <- ggplot(Maximum_RE_DF) + 
  geom_line(
    aes(x = NP + NM, 
    y = RE), 
    color = 'grey',
    linewidth = 0.9) + 
  geom_point(
    aes(x = NP + NM, 
        y = RE, 
        fill = 1-id, 
        group = 1), size = 2.5,
    shape = 21) + 
  geom_line(data = Maximum_Abundance_DF,
    aes(x = NP + NM, 
        y= RE),
        color = '#A6C1A4',
        linewidth= 0.9) + 
  geom_point(data = Maximum_Abundance_DF,
    aes(x = NP + NM, 
        y = RE, 
        fill = 1-id), 
        shape = 22,
        color = 'black',
        size = 2.3) + 
 scale_color_viridis(name = "Primary mortality", limit = c(0,1)) + 
 scale_fill_viridis(name = "Primary mortality", limit = c(0,1)) + 
 xlab(expression("Total vector abundance " * "(" * N[P] *"+" * N[M] * ")")) + 
 ylab(expression(R[E])) + 
 theme_classic() + 
 theme(axis.text = element_text(size = 14, color = 'black'),
       axis.title = element_text(size = 14.5)); GG_Totalvectorabundance_RE

ggsave(here("Figures_Process", "DEC_10","GG_Totalvectorabundance_RE.pdf"),
  width = 6, height = 5, units = 'in' )

GG_Totalvectorabundance_RE <- ggplot(Maximum_RE_DF) + 
  geom_line(
          aes(x = NP + NM, 
              y = RE), 
          color = 'grey',
          linewidth = 0.9) + 
  geom_point(
          aes(x = NP + NM, 
              y = RE, 
              fill = 1-id, 
              group = 1), size = 2.5,
          shape = 21) + 
  geom_line(data = Maximum_Abundance_DF,
            aes(x = NP + NM, 
                y= RE),
            color = '#A6C1A4',
            linewidth= 0.9) + 
  geom_point(data = Maximum_Abundance_DF,
             aes(x = NP + NM, 
                 y = RE, 
                 fill = 1-id), 
             shape = 22,
             color = 'black',
             size = 2.3) + 
  scale_color_viridis(name = "Primary mortality", limit = c(0,1)) + 
  scale_fill_viridis(name = "Primary mortality", limit = c(0,1)) + 
  xlab(expression("Total vector abundance " * "(" * N[P] *"+" * N[M] * ")")) + 
  ylab(expression(R[E])) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 14.5)); GG_Totalvectorabundance_  

   
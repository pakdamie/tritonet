###Assuming that 5% of the total population is infected...

NP = seq(1,2500,10) # Total primary vectors
NM = seq(1,2500,10) # Total secondary vectors 
abund_expand_DF <- expand.grid(NP = NP, NM = NM)

#We are interested in situations where the secondary vector is (1) significantly worse,
# (2) slightly worse (standard), (3) the same, and (4) better. 

params_interest <- c("worse_m","standard", "no_diff", "better_m")

###Calculate the R0 for the different abundances depending on the parameter.
df_expand_RE <- lapply(params_interest, function(x) {
        tmp <- Calculate_Human_Reff_Expanded_static(abund_expand_DF, get_parameters(x))
        tmp$id = x
        
        return(tmp)
        })

df_expand_RE <- do.call(rbind, df_expand_RE)

ggplot(df_expand_RE, aes(x= NP,y = NM, fill = MtoH/RE)) + 
        geom_tile() + facet_wrap(~id) + 
        scale_fill_viridis(option = 'turbo')




saveRDS(df_expand_RE, file = here("Output","df_expand_RE.rds"))

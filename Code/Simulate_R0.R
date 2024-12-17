###Assuming that 5% of the total population is infected...

NP = seq(50,5000,50) # Total primary vectors
NM = seq(50,5000,50) # Total secondary vectors 
abund_expand_DF <- expand.grid(NP = NP, NM = NM)
df_expand_RE <- Calculate_Human_Reff_Expanded(abund_expand_DF, param_standard)

saveRDS(df_expand_RE, file = here("Output","df_expand_RE.rds"))

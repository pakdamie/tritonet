###Assuming that 5% of the total population is infected...

NP = seq(50,5000,100)
NM = seq(50,5000,100)

PI = NP * 0.05
MI = NM * 0.05 

PS = NP - PI
MS = NM - MI

HS = 1000


abund_expand_DF <- expand.grid(NP = NP, NM = NM)

df_expand_RE <- Calculate_Human_Reff_Expanded(abund_expand_DF, param_standard)


heatmap_RE_GG<-
 ggplot(df_expand_RE, 
 aes(x = NP, y= NM , fill = RE)) + 
 scale_x_continuous(expand = c(0,0)) + 
 scale_y_continuous(expand = c(0,0)) + 
 geom_tile() + 
 xlab(expression("Abundance of primary vectors " * "(" * N[P] * ")")) +
 ylab(expression("Abundance of secondary vectors " * "(" * N[M] * ")")) +
 geom_abline(b = 1, color = '#301934') + 
 scale_fill_viridis(option = 'rocket', 
                    name =   expression(R[E])) + 
 theme(axis.text = element_text(size = 12.5, color = 'black'),
       axis.title = element_text(size = 14))


GG_PtoMabundance_RE 
                    
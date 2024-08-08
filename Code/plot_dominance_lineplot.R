plot_dominance_lineplot <- function(calculated_PSV){
        ggplot(calculated_PSV , aes(x = time, y= value, color = patch_num,
                                    group = patch_num))+
                geom_line() + 
                geom_hline(yintercept = 1, color = 'red')+
                scale_color_viridis()+theme_bw()+
                ylab("Species dominance")+
                xlab("Time")+
                theme(legend.position = "none")
        
}
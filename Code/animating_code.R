



### ALL SUBSEQUENT CODE AFTER IS PLACED SOMEWHERE ELSE.
All_individuals = rowSums(results[,(2:ncol(results))])


PSV_df <- calculate_PSV_ratio(results)

df_network <- ggnetwork::fortify(   g9 )
nodes_coord <-  distinct(df_network,x,y,name)

nodes_coord$name <- as.numeric(nodes_coord$name) 

full_PI <- left_join(nodes_coord, PSV_df, by=c("name" = "patch_num"))

a_gif <- ggplot(data = df_network)+
        geom_segment(aes( x=x, xend = xend, y = y, yend = yend))+
        geom_point(data=full_PI,aes(x=x,y=y,color =value),size = 10)+
        scale_colour_gradient2(
                name = "Primary/Secondary",
                low = ("blue"),
                mid = "white",
                high = ("red"),
                midpoint = 1)+
        labs(title = 'Values at {(as.integer(frame_time))}')+
        transition_time(time) + theme_dark()

animate(a_gif, height = 800, width =800)
1anim_save("Gapminder_example.gif")
# X-Y coords
simulated_xy_coord <- simulate_xy_coordinates(24601, 20)
coordinates <- data.frame(simulated_xy_coord [[1]])

ggplot(coordinates, aes(x = x_coord, y = y_coord)) +
  geom_point(size = 3.5, color = "#238A8D") +
  xlab("X") +
  ylab("Y") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 17),
    axis.text = element_blank(),
    axis.ticks = element_blank())

ggsave(here("Figures_Process/Schematic_Figure/xy_coord.pdf"),
       width = 4, height = 3, units = 'in')


# Big components
biggest_component_df <- retrieve_biggest_component(simulated_xy_coord)

ggplot(biggest_component_df,
       aes(x = x, y = y, 
           xend = xend, 
           yend = yend)) +
        geom_edges() + 
        geom_nodes(size = 5.5, 
                   color = '#238A8D') + 
        theme_void()
        
ggsave(here("Figures_Process/Schematic_Figure/network_component.pdf"),
       width = 7, height = 7, units = 'in')
        
